#include "initialization/batch_imu_initializer.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <stdexcept>
using namespace initialization;
void require(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
M3D attitude(double t,bool excited) {
    if(!excited)return M3D::Identity();
    return (Eigen::AngleAxisd(.15*std::sin(.9*t),V3D::UnitZ())*
            Eigen::AngleAxisd(.3*std::sin(1.1*t),V3D::UnitY())*
            Eigen::AngleAxisd(.4*std::sin(.7*t),V3D::UnitX())).toRotationMatrix();
}
std::vector<std::vector<double>> csv(const std::filesystem::path& path) {
    std::ifstream in(path);std::string line;std::getline(in,line);std::vector<std::vector<double>> out;
    while(std::getline(in,line)) {
        std::stringstream row(line);std::string field;std::vector<double> values;
        while(std::getline(row,field,','))values.push_back(std::stod(field));
        if(!values.empty())out.push_back(values);
    }
    return out;
}
int main(int argc,char** argv) {
    if(argc==3) {
        const std::filesystem::path input=argv[1],output=argv[2];
        std::vector<ImuSample::ConstPtr> imu;
        for(const auto& a:csv(input/"imu.csv")) {
            auto s=std::make_shared<ImuSample>();s->timestamp=a[0];
            s->angular_velocity=V3D(a[1],a[2],a[3]);s->acceleration=V3D(a[4],a[5],a[6]);imu.push_back(s);
        }
        std::map<int,std::vector<BatchLidarPose>> segments;
        for(const auto& a:csv(input/"lidar.csv")) {
            const int id=a.size()>11 ? int(a[11]):0;
            segments[id].push_back({a[0],V3D(a[1],a[2],a[3]),Eigen::Quaterniond(a[7],a[4],a[5],a[6]).toRotationMatrix()});
        }
        std::vector<BatchLidarPose> best,poses;
        for(const auto& pair:segments) {
            const auto& block=pair.second;
            if(best.empty() || block.back().time-block.front().time>best.back().time-best.front().time)best=block;
        }
        for(size_t i=0;i<best.size();i+=5)if(best[i].time-best[0].time<=10.)poses.push_back(best[i]);
        const auto result=SolveBatchImuInitialization(poses,imu,output);
        std::cout<<result.report<<'\n';return 0;
    }
    const V3D bg(.003,-.002,.001),ba(.12,-.08,.06),g(.18,-.24,-std::sqrt(9.81*9.81-.18*.18-.24*.24));
    for(bool excited:{false,true}) {
        std::vector<ImuSample::ConstPtr> imu;std::vector<BatchLidarPose> poses;
        for(int k=0;k<=1200;++k) {
            const double t=k*.01,h=1e-5;const M3D R=attitude(t,excited);
            const Eigen::AngleAxisd omega(attitude(t-h,excited).transpose()*attitude(t+h,excited));
            auto s=std::make_shared<ImuSample>();s->timestamp=t;
            s->angular_velocity=omega.axis()*omega.angle()/(2*h)+bg;
            const V3D acceleration(.15*std::sin(.8*t),.1*std::cos(.6*t),.03*std::sin(t));
            s->acceleration=R.transpose()*(acceleration-g)+ba;imu.push_back(s);
            if(k%50==0) {
                const V3D p(.7*t+.15/.8*(t-std::sin(.8*t)/.8),
                             .1/(.6*.6)*(1-std::cos(.6*t)),.03*(t-std::sin(t)));
                poses.push_back({t,p,R});
            }
        }
        const auto result=SolveBatchImuInitialization(poses,imu);
        require(std::abs(result.gravity.norm()-9.81)<1e-10,"gravity norm constraint violated");
        if(!excited) {
            // Soft pose estimates can move slightly under the gyro prior, so
            // the local Jacobian need not contain an exact floating-point null.
            require(result.report["data_only_worst_conditional_std_m_s2"].asDouble()>10.,
                    "constant-attitude ambiguity hidden by priors");
            std::cout<<"PASS constant-attitude data-only degeneracy\n";
        } else {
            require((result.gyro_bias-bg).norm()<.0007,"gyro bias recovery failed");
            require((result.acceleration_bias-ba).norm()<.05,"acceleration bias MAP recovery failed");
            require((result.gravity-g).norm()<.05,"gravity MAP recovery failed");
            require((result.initial_velocity-V3D(.7,0,0)).norm()<.03,"initial velocity recovery failed");
            const std::vector<M3D> frames{
                (Eigen::AngleAxisd(.2,V3D::UnitZ())*Eigen::AngleAxisd(-.1,V3D::UnitY())*Eigen::AngleAxisd(.15,V3D::UnitX())).toRotationMatrix(),
                Eigen::AngleAxisd(3.141592653589793,V3D::UnitX()).toRotationMatrix()};
            for(const auto& Q:frames) {
                auto rotated=poses;for(auto& pose:rotated){pose.position=Q*pose.position;pose.rotation=Q*pose.rotation;}
                const auto other=SolveBatchImuInitialization(rotated,imu);
                require((other.acceleration_bias-result.acceleration_bias).norm()<.0005,"bias depends on world frame");
                require((other.gravity-Q*result.gravity).norm()<.0005,"gravity world-frame transform incorrect");
            }
            std::cout<<"PASS excited bias/gravity/velocity recovery and world-frame invariance\n";
        }
        auto invalid=poses;invalid[1].time=invalid[0].time;
        bool rejected=false;try{SolveBatchImuInitialization(invalid,imu);}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"unordered batch observations accepted");
        if(!excited) {
            auto padded=imu;padded.erase(padded.begin(),padded.begin()+2);
            const auto held=SolveBatchImuInitialization(poses,padded);
            require(std::abs(held.report["leading_imu_hold_s"].asDouble()-.02)<1e-10,"short leading IMU hold not recorded");
            padded.erase(padded.begin(),padded.begin()+20);
            rejected=false;try{SolveBatchImuInitialization(poses,padded);}catch(const std::runtime_error&){rejected=true;}
            require(rejected,"large IMU coverage gap accepted");
        }
    }
    std::cout<<"PASS invalid batch observation rejection\n";
}
