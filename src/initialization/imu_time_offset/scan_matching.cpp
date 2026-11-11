#include "initialization/imu_time_offset/scan_matching.hpp"
#include <pcl/kdtree/kdtree_flann.h>
#include <Eigen/Eigenvalues>
#include <numeric>
#include <stdexcept>

namespace initialization { namespace imu_time_offset {
namespace {
using Policy = livo::ImuTimeOffsetEstimator;
using Matrix6 = Eigen::Matrix<double,6,6>;
using Vector6 = Eigen::Matrix<double,6,1>;
struct Reading { V3D gyro, acceleration; };
Reading ReadImu(const std::vector<ImuSample::ConstPtr>& imu, double t) {
    if (imu.size() < 2 || t < imu.front()->timestamp || t > imu.back()->timestamp)
        throw std::runtime_error("IMU offset deskew would extrapolate");
    auto it = std::upper_bound(imu.begin(),imu.end(),t,
        [](double time,const auto& sample) { return time < sample->timestamp; });
    if (it == imu.end()) return {imu.back()->angular_velocity,imu.back()->acceleration};
    const auto& b = **it; const auto& a = **(it-1);
    const double u = (t-a.timestamp)/(b.timestamp-a.timestamp);
    return {(1-u)*a.angular_velocity+u*b.angular_velocity,
            (1-u)*a.acceleration+u*b.acceleration};
}
}

bool HasImuCoverage(const RotationCloud& input,
                    const std::vector<ImuSample::ConstPtr>& imu, double padding) {
    if (imu.size() < 2 || input.point_times.empty() ||
        input.point_times.size() != input.cloud->size()) return false;
    const auto bounds = std::minmax_element(input.point_times.begin(),input.point_times.end());
    const double begin = std::min(input.time,*bounds.first)-padding;
    const double end = std::max(input.time,*bounds.second)+padding;
    if (begin < imu.front()->timestamp || end > imu.back()->timestamp) return false;
    auto it = std::upper_bound(imu.begin(),imu.end(),begin,
        [](double time,const auto& sample) { return time < sample->timestamp; });
    for (; it != imu.end(); ++it) {
        const auto& a = **(it-1); const auto& b = **it;
        if (b.timestamp <= a.timestamp || b.timestamp-a.timestamp > .025 ||
            !a.acceleration.allFinite() || !a.angular_velocity.allFinite() ||
            !b.acceleration.allFinite() || !b.angular_velocity.allFinite()) return false;
        // Include the sample bracketing the right edge.
        if (b.timestamp >= end) break;
    }
    return true;
}

DeskewState QueryDeskewState(const RotationCloud& input,
                             const std::vector<ImuSample::ConstPtr>& imu,
                             const LidarTwist& twist, double acceleration_scale) {
    DeskewState state;
    state.velocity = twist.velocity;
    state.acceleration_scale = acceleration_scale;
    const auto it = std::lower_bound(imu.begin(),imu.end(),input.time,
        [](const auto& sample,double time) { return sample->timestamp < time; });
    if (it != imu.end() && (*it)->orientation_valid &&
        (*it)->orientation.coeffs().allFinite() && (*it)->orientation.norm() > 1e-6)
        state.gravity = (*it)->orientation.normalized().conjugate()*V3D(0,0,-9.81);
    else {
        // Same configured-clock gravity seed for ALL candidates. This is not
        // refitted to make any individual candidate look better.
        const V3D specific_force = ReadImu(imu,input.time).acceleration*acceleration_scale;
        const V3D seed = twist.acceleration-specific_force;
        if (seed.norm() < 1e-6) throw std::runtime_error("invalid IMU offset gravity seed");
        state.gravity = 9.81*seed.normalized();
    }
    return state;
}

RotationCloud DeskewWithImu(const RotationCloud& input,
                            const std::vector<ImuSample::ConstPtr>& imu,
                            const DeskewState& state, double shift) {
    if (!HasImuCoverage(input,imu,std::abs(shift)))
        throw std::runtime_error("IMU offset deskew lacks common timestamp coverage");
    RotationCloud out; out.time = input.time; out.point_times = input.point_times;
    *out.cloud = *input.cloud;
    std::vector<size_t> order(input.cloud->size());
    std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](size_t a,size_t b) {
        return input.point_times[a] < input.point_times[b];
    });
    const auto split = std::lower_bound(order.begin(),order.end(),input.time,
        [&](size_t index,double t) { return input.point_times[index] < t; });
    auto integrate = [&](auto first,auto last) {
        M3D rotation = M3D::Identity(); V3D position = V3D::Zero(), velocity = state.velocity;
        double time = input.time;
        for (auto it = first; it != last; ++it) {
            const size_t index = *it; const double target = input.point_times[index];
            while (std::abs(target-time) > 1e-12) {
                const double dt = std::clamp(target-time,-.001,.001);
                // Positive shift adds to tags: effective(t) reads raw(t-shift).
                const auto a = ReadImu(imu,time-shift), b = ReadImu(imu,time+dt-shift);
                const V3D w = .5*(a.gyro+b.gyro)-state.gyro_bias;
                const M3D half = rotation*Policy::Exp(w*(.5*dt));
                const V3D acceleration = half*(.5*(a.acceleration+b.acceleration)*
                    state.acceleration_scale-state.acceleration_bias)+state.gravity;
                position += velocity*dt+.5*acceleration*dt*dt;
                velocity += acceleration*dt;
                rotation = rotation*Policy::Exp(w*dt); time += dt;
            }
            const V3D corrected = rotation*(*input.cloud)[index].getVector3fMap().cast<double>()+position;
            auto& point = (*out.cloud)[index];
            point.x = corrected.x(); point.y = corrected.y(); point.z = corrected.z();
        }
    };
    integrate(split,order.end());
    integrate(std::make_reverse_iterator(split),order.rend());
    return out;
}

struct FrozenPlaneMap::Impl {
    PointCloudXYZI::Ptr cloud;
    pcl::KdTreeFLANN<PointType> tree;
    std::vector<V3D> normals, centers;
    explicit Impl(const RotationCloud& target) : cloud(new PointCloudXYZI(*target.cloud)),
        normals(cloud->size(),V3D::Zero()), centers(cloud->size(),V3D::Zero()) {
        if (cloud->size() < 12) return;
        tree.setInputCloud(cloud);
        for (size_t i = 0; i < cloud->size(); ++i) {
            std::vector<int> ids(12);std::vector<float> distances(12);
            if (tree.nearestKSearch((*cloud)[i],12,ids,distances) != 12 || distances.back() > 2.25) continue;
            V3D mean = V3D::Zero();
            for (int id : ids) mean += (*cloud)[id].getVector3fMap().cast<double>();
            mean /= 12.; M3D covariance = M3D::Zero();
            for (int id : ids) {
                const V3D d = (*cloud)[id].getVector3fMap().cast<double>()-mean;
                covariance += d*d.transpose()/12.;
            }
            Eigen::SelfAdjointEigenSolver<M3D> es(covariance);
            if (es.eigenvalues()[1] < .005 || es.eigenvalues()[0] > .01 ||
                es.eigenvalues()[0] > .15*es.eigenvalues()[1]) continue;
            normals[i] = es.eigenvectors().col(0); centers[i] = mean;
        }
    }
    bool Association(const V3D& point, V3D& normal, double& residual) const {
        if (cloud->size() < 12) return false;
        PointType search{}; search.x=point.x();search.y=point.y();search.z=point.z();
        std::vector<int> ids(1);std::vector<float> distances(1);
        if (tree.nearestKSearch(search,1,ids,distances) != 1 || distances[0] > .64 ||
            normals[ids[0]].norm() < .5) return false;
        normal = normals[ids[0]];
        residual = normal.dot(point-centers[ids[0]]);
        return std::isfinite(residual);
    }
};
FrozenPlaneMap::FrozenPlaneMap(const RotationCloud& target) : impl_(new Impl(target)) {}
FrozenPlaneMap::~FrozenPlaneMap() = default;
FrozenPlaneMap::FrozenPlaneMap(FrozenPlaneMap&&) noexcept = default;
FrozenPlaneMap& FrozenPlaneMap::operator=(FrozenPlaneMap&&) noexcept = default;

livo::OffsetMatchScore FrozenPlaneMap::Evaluate(const RotationCloud& source) const {
    livo::OffsetMatchScore score; score.scans = 1; score.points = source.cloud->size();
    const double cap2 = Policy::kResidualCap*Policy::kResidualCap;
    score.capped_squared_error = cap2*score.points;
    M3D R = M3D::Identity(); V3D t = V3D::Zero();
    // Same optimizer and identity prior for every candidate; no covariance or
    // candidate-dependent gate can increase the denominator or loosen scoring.
    for (int iteration = 0; iteration < 8; ++iteration) {
        Matrix6 H = Matrix6::Zero(); Vector6 g = Vector6::Zero(); size_t count = 0;
        for (const auto& point : *source.cloud) {
            const V3D rotated = R*point.getVector3fMap().cast<double>(), p = rotated+t;
            V3D n; double e;
            if (!impl_->Association(p,n,e) || std::abs(e) > Policy::kResidualCap) continue;
            const double w = std::min(1.,.05/std::max(std::abs(e),1e-12));
            Vector6 J; J << rotated.cross(n),n;
            H += w*J*J.transpose(); g += w*J*e; ++count;
        }
        if (count < 100) return score;
        const Vector6 step = -(H+Matrix6::Identity()*1e-6).ldlt().solve(g);
        if (!step.allFinite() || step.head<3>().norm() > .0873 || step.tail<3>().norm() > .75) return score;
        R = Policy::Exp(step.head<3>())*R; t += step.tail<3>();
        if (step.head<3>().norm() < 1e-5 && step.tail<3>().norm() < 1e-4) break;
    }
    // Fresh associations at FINAL pose. Unmatched points remain in both the
    // fixed ratio denominator and the capped-residual score.
    score.valid_scans = 1; score.capped_squared_error = 0;
    for (const auto& point : *source.cloud) {
        const V3D p = R*point.getVector3fMap().cast<double>()+t;
        V3D n; double e;
        if (!impl_->Association(p,n,e)) { score.capped_squared_error += cap2; continue; }
        score.capped_squared_error += std::min(cap2,e*e);
        if (std::abs(e) <= Policy::kMatchDistance) {
            ++score.matches; score.squared_error += e*e;
        }
    }
    return score;
}
} }  // namespace initialization::imu_time_offset
