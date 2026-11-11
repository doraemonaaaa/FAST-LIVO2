#include "initialization/imu_time_offset/lidar_rotation.hpp"
#include "initialization/imu_time_offset/estimator.hpp"
#include <pcl/kdtree/kdtree_flann.h>
#include <Eigen/Eigenvalues>
#include <map>

namespace initialization { namespace imu_time_offset {
using Matrix6 = Eigen::Matrix<double,6,6>;
using Vector6 = Eigen::Matrix<double,6,1>;

// Select an acquired point per voxel, retaining its own position and time.
// A centroid of points acquired at different times is not a timed observation.
RotationPhases RotationClouds(const LivoScan& scan, const Eigen::Matrix4d& T) {
    struct Cell { V3D point; double time, distance; };
    std::array<std::map<std::array<int,3>,Cell>,kPhases> grids;
    for (const auto& point : scan.points) {
        const int index = int(std::floor((point.time-scan.begin)/kSliceSeconds));
        if (index < 0 || index >= kPhases) continue;
        const V3D p = T.block<3,3>(0,0)*V3D(point.point.x,point.point.y,point.point.z)
                     +T.block<3,1>(0,3);
        if (p.norm() < 1 || p.norm() > 45) continue;
        const std::array<int,3> key{{int(std::floor(p.x()/.25)),
            int(std::floor(p.y()/.25)),int(std::floor(p.z()/.25))}};
        const V3D center((key[0]+.5)*.25,(key[1]+.5)*.25,(key[2]+.5)*.25);
        const double distance = (p-center).squaredNorm();
        auto it = grids[index].find(key);
        if (it == grids[index].end() || distance < it->second.distance)
            grids[index][key] = {p,point.time,distance};
    }
    RotationPhases result;
    for (int i = 0; i < kPhases; ++i) {
        result[i].time = scan.begin+(i+.5)*kSliceSeconds;
        if (scan.end < scan.begin+(i+1)*kSliceSeconds-.001) continue;
        for (const auto& entry : grids[i]) {
            const auto& cell = entry.second;
            PointType point{};point.x=cell.point.x();point.y=cell.point.y();point.z=cell.point.z();
            result[i].cloud->push_back(point);
            result[i].point_times.push_back(cell.time);
        }
    }
    return result;
}

RotationCloud FullRotationCloud(const RotationPhases& phases) {
    RotationCloud out;
    out.time = .5*(phases.front().time+phases.back().time);
    struct Cell { PointType point; double time, distance; };
    std::map<std::array<int,3>,Cell> grid;
    for (const auto& phase : phases) {
        for (size_t i = 0; i < phase.cloud->size(); ++i) {
            const auto& point = (*phase.cloud)[i];
            const V3D p = point.getVector3fMap().cast<double>();
            const std::array<int,3> key{{int(std::floor(p.x()/.25)),
                int(std::floor(p.y()/.25)),int(std::floor(p.z()/.25))}};
            const V3D center((key[0]+.5)*.25,(key[1]+.5)*.25,(key[2]+.5)*.25);
            const double distance = (p-center).squaredNorm();
            auto it = grid.find(key);
            if (it == grid.end() || distance < it->second.distance)
                grid[key] = {point,phase.point_times[i],distance};
        }
    }
    for (const auto& entry : grid) {
        out.cloud->push_back(entry.second.point);
        out.point_times.push_back(entry.second.time);
    }
    return out;
}

LidarTwist EstimateTwist(const RotationFit& backward, double dt_backward,
                         const RotationFit& forward, double dt_forward) {
    LidarTwist twist;
    // Both sides are required: one-sided interval rates lag/lead the reference
    // by half a scan, precisely the effect a temporal calibration must avoid.
    if (!backward.valid || !forward.valid || dt_backward <= 0 || dt_forward <= 0)
        return twist;
    const V3D wb = livo::ImuTimeOffsetEstimator::Log(backward.rotation)/dt_backward;
    const V3D wf = livo::ImuTimeOffsetEstimator::Log(forward.rotation)/dt_forward;
    const V3D vb = backward.rotation.transpose()*backward.translation/dt_backward;
    const V3D vf = forward.translation/dt_forward;
    const double span = dt_backward+dt_forward;
    twist.angular_velocity = (dt_forward*wb+dt_backward*wf)/span;
    twist.velocity = (dt_forward*vb+dt_backward*vf)/span;
    twist.angular_acceleration = 2*(wf-wb)/span;
    twist.acceleration = 2*(vf-vb)/span;
    twist.valid = twist.angular_velocity.allFinite() && twist.velocity.allFinite() &&
        twist.angular_acceleration.allFinite() && twist.acceleration.allFinite();
    return twist;
}

LidarTwist ShiftTwistReference(const LidarTwist& input, double dt) {
    if (!input.valid) return input;
    const M3D rotation = livo::ImuTimeOffsetEstimator::Exp(
        (input.angular_velocity*dt+.5*input.angular_acceleration*dt*dt).eval());
    LidarTwist out = input;
    out.angular_velocity = rotation.transpose()*(input.angular_velocity+
                                                  input.angular_acceleration*dt);
    out.velocity = rotation.transpose()*(input.velocity+input.acceleration*dt);
    out.angular_acceleration = rotation.transpose()*input.angular_acceleration;
    out.acceleration = rotation.transpose()*input.acceleration;
    return out;
}

RotationCloud DeskewRotationCloud(const RotationCloud& input, const LidarTwist& twist) {
    RotationCloud out;out.time = input.time;
    if (!twist.valid || input.point_times.size() != input.cloud->size()) return out;
    out.point_times = input.point_times;out.cloud->reserve(input.cloud->size());
    for (size_t i = 0; i < input.cloud->size(); ++i) {
        const double dt = input.point_times[i]-input.time;
        const V3D p = (*input.cloud)[i].getVector3fMap().cast<double>();
        const V3D corrected = livo::ImuTimeOffsetEstimator::Exp(
            (twist.angular_velocity*dt+.5*twist.angular_acceleration*dt*dt).eval())*p
            +twist.velocity*dt+.5*twist.acceleration*dt*dt;
        PointType point{};point.x=corrected.x();point.y=corrected.y();point.z=corrected.z();
        out.cloud->push_back(point);
    }
    return out;
}

// LiDAR-only point-to-plane alignment; no gyro prior and no IMU deskew.
RotationFit RegisterRotation(const RotationCloud& target,
                              const RotationCloud& source,
                              const RotationFit& prediction) {
    RotationFit out;
    if (target.cloud->size() < 100 || source.cloud->size() < 100) return out;
    pcl::KdTreeFLANN<PointType> tree;
    tree.setInputCloud(target.cloud);
    std::vector<V3D> normals(target.cloud->size(),V3D::Zero());
    std::vector<V3D> centers(target.cloud->size(),V3D::Zero());
    for (size_t i = 0; i < target.cloud->size(); ++i) {
        std::vector<int> ids(12);std::vector<float> distances(12);
        if (tree.nearestKSearch((*target.cloud)[i],12,ids,distances) != 12 ||
            distances.back() > 2.25) continue;
        V3D mean = V3D::Zero();
        for (int id : ids) mean += (*target.cloud)[id].getVector3fMap().cast<double>();
        mean /= 12.;M3D covariance = M3D::Zero();
        for (int id : ids) {
            const V3D d = (*target.cloud)[id].getVector3fMap().cast<double>()-mean;
            covariance += d*d.transpose()/12.;
        }
        Eigen::SelfAdjointEigenSolver<M3D> es(covariance);
        if (es.eigenvalues()[1] < .005 || es.eigenvalues()[0] > .01 ||
            es.eigenvalues()[0] > .15*es.eigenvalues()[1]) continue;
        normals[i] = es.eigenvectors().col(0);centers[i] = mean;
    }
    M3D R = prediction.valid ? prediction.rotation : M3D::Identity();
    V3D t = prediction.valid ? prediction.translation : V3D::Zero();
    auto linearize = [&](Matrix6& H,Vector6& g,double& squared,size_t& count) {
        H.setZero();g.setZero();squared=0;count=0;
        std::vector<int> ids(1);std::vector<float> distances(1);
        for (const auto& point : *source.cloud) {
            const V3D rotated = R*V3D(point.x,point.y,point.z),p = rotated+t;
            PointType search{};search.x=p.x();search.y=p.y();search.z=p.z();
            if (tree.nearestKSearch(search,1,ids,distances) != 1 ||
                distances[0] > .64 || normals[ids[0]].norm() < .5) continue;
            const V3D& n = normals[ids[0]];
            const double e = n.dot(p-centers[ids[0]]);
            if (std::abs(e) > .3) continue;
            const double w = std::min(1.,.05/std::max(std::abs(e),1e-12));
            Vector6 J;J << rotated.cross(n),n;
            H += w*J*J.transpose();g += w*J*e;squared += e*e;++count;
        }
    };
    Matrix6 H;Vector6 g;double squared;size_t count;
    for (int iteration = 0; iteration < 25; ++iteration) {
        linearize(H,g,squared,count);
        if (count < 100) return out;
        const Vector6 step = -(H+Matrix6::Identity()*1e-9).ldlt().solve(g);
        if (!step.allFinite() || step.head<3>().norm() > .0873 ||
            step.tail<3>().norm() > .75) return out;
        R = livo::ImuTimeOffsetEstimator::Exp(step.head<3>())*R;t += step.tail<3>();
        if (step.head<3>().norm() < 1e-5 && step.tail<3>().norm() < 1e-4) break;
    }
    linearize(H,g,squared,count);
    out.matches = count;out.ratio = double(count)/source.cloud->size();
    out.rms = count ? std::sqrt(squared/count) : 0;
    if (count < 100 || out.ratio < .25 || out.rms > .10) return out;
    const M3D Htt = H.bottomRightCorner<3,3>();
    Eigen::SelfAdjointEigenSolver<M3D> translation_geometry(Htt);
    if (translation_geometry.eigenvalues().minCoeff() < 1e-5) return out;
    const M3D Htr = H.bottomLeftCorner<3,3>();
    const M3D information = (H.topLeftCorner<3,3>()-
        Htr.transpose()*Htt.ldlt().solve(Htr)).eval();
    Eigen::SelfAdjointEigenSolver<M3D> es((.5*(information+information.transpose())).eval());
    if (es.eigenvalues().minCoeff() < 1e-5 ||
        es.eigenvalues().minCoeff()/es.eigenvalues().maxCoeff() < 1e-4) return out;
    const M3D covariance = information.ldlt().solve(M3D::Identity())*
        std::max(.03*.03,out.rms*out.rms);
    if (!covariance.allFinite() || covariance.diagonal().maxCoeff() > .00524*.00524)
        return out;
    out.valid = true;out.rotation = R;out.translation = t;out.covariance = covariance;
    return out;
}
} }  // namespace initialization::imu_time_offset
