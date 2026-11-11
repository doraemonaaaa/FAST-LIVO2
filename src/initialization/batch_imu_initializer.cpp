#include "initialization/batch_imu_initializer.hpp"
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace initialization { namespace {
using Q=Eigen::Quaterniond;
using V2=Eigen::Vector2d;
using M32=Eigen::Matrix<double,3,2>;
constexpr double kGravity=9.81;
Q Exp(const V3D& v) {
    double q[4];ceres::AngleAxisToQuaternion(v.data(),q);return Q(q[0],q[1],q[2],q[3]);
}
V3D Log(const Q& q) {
    const double data[4]={q.w(),q.x(),q.y(),q.z()};V3D out;
    ceres::QuaternionToAngleAxis(data,out.data());return out;
}
struct GravityChart {
    int axis=2,first=0,second=1;
    double sign=-1.;
    V3D coordinates(const double* u) const {
        V3D q=V3D::Zero();q[axis]=sign;q[first]=u[0];q[second]=u[1];return q;
    }
    V3D value(const double* u) const {return kGravity*coordinates(u).normalized();}
};

struct Reading {V3D gyro,acceleration;};
Reading Interpolate(const std::vector<ImuSample::ConstPtr>& imu,double time) {
    auto next=std::lower_bound(imu.begin(),imu.end(),time,[](const auto& s,double t){return s->timestamp<t;});
    if(next==imu.end()) return {imu.back()->angular_velocity,imu.back()->acceleration};
    if(next==imu.begin()) return {(*next)->angular_velocity,(*next)->acceleration};
    const auto& a=**(next-1);const auto& b=**next;
    const double f=(time-a.timestamp)/(b.timestamp-a.timestamp);
    return {(1.-f)*a.angular_velocity+f*b.angular_velocity,(1.-f)*a.acceleration+f*b.acceleration};
}
struct Interval {
    double dt=0;
    Q rotation=Q::Identity();
    V3D dv=V3D::Zero(),dp=V3D::Zero();
    M3D Bv=M3D::Zero(),Bp=M3D::Zero();
    M3D Jr=M3D::Zero(),Jv=M3D::Zero(),Jp=M3D::Zero();
};
Interval Integrate(const std::vector<ImuSample::ConstPtr>& imu,double begin,double end,const V3D& bg) {
    // A scan midpoint may precede the first IMU by part of one scan. Match
    // the reviewed interpolation's endpoint hold, but bound missing coverage.
    if(begin<imu.front()->timestamp-.1 || end>imu.back()->timestamp+.1 || end<=begin)
        throw std::runtime_error("batch calibration has insufficient IMU coverage");
    std::vector<double> times{begin};
    auto next=std::upper_bound(imu.begin(),imu.end(),begin,[](double t,const auto& s){return t<s->timestamp;});
    for(;next!=imu.end() && (*next)->timestamp<end;++next) times.push_back((*next)->timestamp);
    times.push_back(end);Interval out;out.dt=end-begin;
    auto previous=Interpolate(imu,begin);
    for(size_t i=1;i<times.size();++i) {
        const double dt=times[i]-times[i-1];const auto current=Interpolate(imu,times[i]);
        const V3D omega=.5*(previous.gyro+current.gyro)-bg;
        const M3D middle=(out.rotation*Exp(omega*(.5*dt))).toRotationMatrix();
        const V3D acceleration=middle*(.5*(previous.acceleration+current.acceleration));
        out.dp+=out.dv*dt+.5*acceleration*dt*dt;
        out.Bp+=out.Bv*dt+.5*middle*dt*dt;
        out.dv+=acceleration*dt;out.Bv+=middle*dt;
        out.rotation=out.rotation*Exp(omega*dt);previous=current;
    }
    return out;
}
Interval Linearize(const std::vector<ImuSample::ConstPtr>& imu,double begin,double end,const V3D& bg) {
    auto out=Integrate(imu,begin,end,bg);constexpr double epsilon=1e-5;
    for(int j=0;j<3;++j) {
        V3D shifted=bg;shifted[j]+=epsilon;
        const auto other=Integrate(imu,begin,end,shifted);
        out.Jr.col(j)=Log(out.rotation.conjugate()*other.rotation)/epsilon;
        out.Jv.col(j)=(other.dv-out.dv)/epsilon;
        out.Jp.col(j)=(other.dp-out.dp)/epsilon;
    }
    return out;
}
struct GyroResidual {
    Interval interval;Q observed;
    bool operator()(const double* bg,double* r) const {
        const Q predicted=interval.rotation*Exp(interval.Jr*Eigen::Map<const V3D>(bg));
        Eigen::Map<V3D> residual(r);residual=Log(predicted.conjugate()*observed)/.01;return true;
    }
};
struct Prior {
    double sigma;
    bool operator()(const double* p,double* r) const {
        Eigen::Map<V3D> residual(r);residual=Eigen::Map<const V3D>(p)/sigma;return true;
    }
};
struct State {V3D correction=V3D::Zero(),velocity=V3D::Zero(),position=V3D::Zero();};
struct JointResidual {
    Interval interval;Q reference_i,reference_j,measured_rotation;
    V3D measured_translation,bg_reference;
    GravityChart chart;
    bool operator()(const double* theta_i,const double* theta_j,const double* vi,const double* vj,
                    const double* pi,const double* pj,const double* bg,const double* ba,const double* u,
                    double* residuals) const {
        const Q ri=Exp(Eigen::Map<const V3D>(theta_i))*reference_i;
        const Q rj=Exp(Eigen::Map<const V3D>(theta_j))*reference_j;
        const V3D dbg=Eigen::Map<const V3D>(bg)-bg_reference;
        const V3D bias=Eigen::Map<const V3D>(ba),g=chart.value(u);
        const Q dr=interval.rotation*Exp(interval.Jr*dbg);
        const V3D dv=interval.dv+interval.Jv*dbg-interval.Bv*bias;
        const V3D dp=interval.dp+interval.Jp*dbg-interval.Bp*bias;
        const double dt=interval.dt;
        const V3D velocity_i=Eigen::Map<const V3D>(vi),velocity_j=Eigen::Map<const V3D>(vj);
        const V3D displacement=Eigen::Map<const V3D>(pj)-Eigen::Map<const V3D>(pi);
        Eigen::Map<Eigen::Matrix<double,15,1>> r(residuals);
        r.segment<3>(0)=Log(dr.conjugate()*(ri.conjugate()*rj))/(.01*std::sqrt(dt));
        r.segment<3>(3)=(ri.conjugate()*(velocity_j-velocity_i-g*dt)-dv)/(.05*std::sqrt(dt));
        r.segment<3>(6)=(ri.conjugate()*(displacement-velocity_i*dt-.5*g*dt*dt)-dp)/(.03*dt);
        r.segment<3>(9)=Log(measured_rotation.conjugate()*(ri.conjugate()*rj))/(.005*std::sqrt(dt/.1));
        r.segment<3>(12)=(ri.conjugate()*displacement-measured_translation)/(.03*std::sqrt(dt/.1));
        return r.allFinite();
    }
};
ceres::Solver::Options Options() {
    ceres::Solver::Options options;options.linear_solver_type=ceres::DENSE_QR;
    options.max_num_iterations=180;options.num_threads=1;
    options.function_tolerance=1e-10;options.gradient_tolerance=1e-10;options.parameter_tolerance=1e-10;
    return options;
}
void Append(Json::Value& out,const V3D& v) {for(int i=0;i<3;++i)out.append(v[i]);}
}  // namespace

BatchImuInitialization SolveBatchImuInitialization(const std::vector<BatchLidarPose>& poses,
    const std::vector<ImuSample::ConstPtr>& imu,const std::filesystem::path& diagnostics) {
    if(poses.size()<6 || imu.size()<2) throw std::runtime_error("insufficient batch calibration data");
    for(size_t i=0;i<imu.size();++i)
        if(!imu[i] || !std::isfinite(imu[i]->timestamp) || !imu[i]->angular_velocity.allFinite() ||
           !imu[i]->acceleration.allFinite() || (i && imu[i]->timestamp<=imu[i-1]->timestamp))
            throw std::runtime_error("invalid batch calibration IMU samples");
    for(size_t i=0;i<poses.size();++i)
        if(!std::isfinite(poses[i].time) || !poses[i].position.allFinite() || !poses[i].rotation.allFinite() ||
           (i && poses[i].time<=poses[i-1].time)) throw std::runtime_error("invalid batch LiDAR trajectory");
    BatchImuInitialization out;out.calibration_begin=poses.front().time;out.calibration_end=poses.back().time;
    std::vector<Interval> intervals;
    for(size_t i=0;i+1<poses.size();++i) intervals.push_back(Linearize(imu,poses[i].time,poses[i+1].time,V3D::Zero()));
    V3D bg=V3D::Zero();ceres::Problem gyro_problem;
    for(size_t i=0;i<intervals.size();++i) {
        auto* cost=new ceres::NumericDiffCostFunction<GyroResidual,ceres::CENTRAL,3,3>(
            new GyroResidual{intervals[i],Q(poses[i].rotation.transpose()*poses[i+1].rotation)});
        gyro_problem.AddResidualBlock(cost,nullptr,bg.data());
    }
    gyro_problem.AddResidualBlock(new ceres::NumericDiffCostFunction<Prior,ceres::CENTRAL,3,3>(new Prior{.01}),nullptr,bg.data());
    ceres::Solver::Summary gyro_summary;ceres::Solve(Options(),&gyro_problem,&gyro_summary);
    if(!gyro_summary.IsSolutionUsable() || !bg.allFinite()) throw std::runtime_error("batch gyro bias fit failed");
    const V3D bg_reference=bg;
    for(size_t i=0;i<intervals.size();++i) intervals[i]=Linearize(imu,poses[i].time,poses[i+1].time,bg_reference);
    std::vector<State> states(poses.size());
    for(size_t i=0;i<poses.size();++i) {
        states[i].position=poses[i].position-poses.front().position;
        const size_t a=i ? i-1:i,b=std::min(i+1,poses.size()-1);
        states[i].velocity=(poses[b].position-poses[a].position)/(poses[b].time-poses[a].time);
    }
    V3D ba=V3D::Zero();V2 gravity_chart=V2::Zero();
    // Acceleration supplies an initial guess only, never an attitude factor.
    const V3D guess=-poses.front().rotation*Interpolate(imu,poses.front().time).acceleration;
    GravityChart chart;
    guess.cwiseAbs().maxCoeff(&chart.axis);
    chart.sign=guess[chart.axis]>=0 ? 1.:-1.;
    chart.first=(chart.axis+1)%3;chart.second=(chart.axis+2)%3;
    if(std::abs(guess[chart.axis])>1e-6)
        gravity_chart=V2(guess[chart.first],guess[chart.second])/std::abs(guess[chart.axis]);
    ceres::Problem problem;std::vector<ceres::ResidualBlockId> data_blocks;
    for(size_t i=0;i<intervals.size();++i) {
        auto* f=new JointResidual{intervals[i],Q(poses[i].rotation),Q(poses[i+1].rotation),
            Q(poses[i].rotation.transpose()*poses[i+1].rotation),
            poses[i].rotation.transpose()*(poses[i+1].position-poses[i].position),bg_reference,chart};
        auto* cost=new ceres::NumericDiffCostFunction<JointResidual,ceres::CENTRAL,15,3,3,3,3,3,3,3,3,2>(f);
        data_blocks.push_back(problem.AddResidualBlock(cost,nullptr,
            states[i].correction.data(),states[i+1].correction.data(),states[i].velocity.data(),states[i+1].velocity.data(),
            states[i].position.data(),states[i+1].position.data(),bg.data(),ba.data(),gravity_chart.data()));
    }
    problem.SetParameterBlockConstant(states[0].correction.data());problem.SetParameterBlockConstant(states[0].position.data());
    problem.AddResidualBlock(new ceres::NumericDiffCostFunction<Prior,ceres::CENTRAL,3,3>(new Prior{.01}),nullptr,bg.data());
    problem.AddResidualBlock(new ceres::NumericDiffCostFunction<Prior,ceres::CENTRAL,3,3>(new Prior{.3}),nullptr,ba.data());
    ceres::Solver::Summary summary;ceres::Solve(Options(),&problem,&summary);
    if(!summary.IsSolutionUsable() || !bg.allFinite() || !ba.allFinite() || !gravity_chart.allFinite())
        throw std::runtime_error("joint LiDAR-IMU batch initialization failed: "+summary.BriefReport());
    out.gyro_bias=bg;out.acceleration_bias=ba;out.gravity=chart.value(gravity_chart.data());out.initial_velocity=states.front().velocity;
    auto& report=out.report;
    report["method"]="joint_lidar_imu_batch_bias_initialization";
    report["calibration_begin_s"]=out.calibration_begin;report["calibration_end_s"]=out.calibration_end;
    report["leading_imu_hold_s"]=std::max(0.,imu.front()->timestamp-out.calibration_begin);
    report["trailing_imu_hold_s"]=std::max(0.,out.calibration_end-imu.back()->timestamp);
    report["keyframes"]=Json::UInt64(poses.size());report["gyro_bias_prior_sigma_rad_s"]=.01;
    report["acceleration_bias_prior_sigma_m_s2"]=.3;report["gravity_norm_m_s2"]=out.gravity.norm();
    report["solver"]=summary.BriefReport();report["initial_cost"]=summary.initial_cost;report["final_cost"]=summary.final_cost;
    report["recorded_lio_used"]=false;report["bias_continues_updating"]=true;
    report["handoff"]="bias means only; original ESIKF prior; refit early velocity/gravity before replay";
    Append(report["gyro_bias"],bg);Append(report["acceleration_bias"],ba);
    Append(report["gravity_in_calibration_reference"],out.gravity);Append(report["initial_velocity_in_calibration_reference"],out.initial_velocity);
    // Marginalize the pose/velocity/gyro variables using DATA factors only.
    // These conditional diagnostics do not claim physical bias observability.
    ceres::Problem::EvaluateOptions evaluate;evaluate.residual_blocks=data_blocks;
    for(size_t i=0;i<states.size();++i) {
        if(i) evaluate.parameter_blocks.push_back(states[i].correction.data());
        evaluate.parameter_blocks.push_back(states[i].velocity.data());
        if(i) evaluate.parameter_blocks.push_back(states[i].position.data());
    }
    evaluate.parameter_blocks.push_back(bg.data());evaluate.parameter_blocks.push_back(ba.data());evaluate.parameter_blocks.push_back(gravity_chart.data());
    ceres::CRSMatrix sparse;
    if(problem.Evaluate(evaluate,nullptr,nullptr,nullptr,&sparse)) {
        Eigen::MatrixXd J=Eigen::MatrixXd::Zero(sparse.num_rows,sparse.num_cols);
        for(int row=0;row<sparse.num_rows;++row)
            for(int j=sparse.rows[row];j<sparse.rows[row+1];++j) J(row,sparse.cols[j])=sparse.values[j];
        const Eigen::MatrixXd nuisance=J.leftCols(J.cols()-5);
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(nuisance);
        const Eigen::MatrixXd Qn=qr.householderQ()*Eigen::MatrixXd::Identity(J.rows(),qr.rank());
        const V3D q=chart.coordinates(gravity_chart.data()),unit=q.normalized();
        const M3D derivative=kGravity/q.norm()*(M3D::Identity()-unit*unit.transpose());
        M32 Jg;Jg.col(0)=derivative.col(chart.first);Jg.col(1)=derivative.col(chart.second);
        M32 tangent;tangent.col(0)=unit.unitOrthogonal();tangent.col(1)=unit.cross(tangent.col(0));
        Eigen::MatrixXd calibration(J.rows(),5);calibration.leftCols<3>()=J.middleCols(J.cols()-5,3);
        calibration.rightCols<2>()=J.rightCols<2>()*Jg.colPivHouseholderQr().solve(tangent);
        const Eigen::MatrixXd projected=calibration-Qn*(Qn.transpose()*calibration);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(projected,Eigen::ComputeThinV);
        for(int i=0;i<5;++i) report["data_only_singular_values"].append(svd.singularValues()[i]);
        if(svd.singularValues()[4]>1e-12) report["data_only_worst_conditional_std_m_s2"]=1./svd.singularValues()[4];
        report["observability_caveat"]="local conditional uncertainty under assumed LiDAR/IMU noise; priors excluded; weak bias/gravity directions may remain";
    }
    if(!diagnostics.empty()) {
        std::filesystem::create_directories(diagnostics);
        std::ofstream file(diagnostics/"batch_lidar_poses.csv");file<<std::setprecision(17)<<"time_s,x,y,z,qx,qy,qz,qw\n";
        for(const auto& p:poses) {
            const Q q(p.rotation);file<<p.time;
            for(int j=0;j<3;++j)file<<','<<p.position[j];
            for(int j=0;j<4;++j)file<<','<<q.coeffs()[j];
            file<<'\n';
        }
        std::ofstream(diagnostics/"batch_imu_initialization.json")<<report;
    }
    return out;
}
}
