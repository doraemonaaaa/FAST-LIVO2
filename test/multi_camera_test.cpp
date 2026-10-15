#include "vio.h"
#include "state_numerics.h"
#include "camera_sync.h"
#include <iostream>
#include <stdexcept>

static void require(bool result, const char *message)
{
  if (!result) throw std::runtime_error(message);
}

int main()
{
  try
  {
    MeasureGroup bundle;
    VisualPoint shared_point(V3D::Zero(),4);
    shared_point.referenceConverged(0)=true;
    require(!shared_point.referenceConverged(1),"one camera cannot freeze another camera's reference lifecycle");
    bundle.images={{0,1,cv::Mat(2,2,CV_8UC1,cv::Scalar(7)),1000000000},
                   {1,1,cv::Mat(2,2,CV_8UC1,cv::Scalar(9)),1000000000}};
    bundle.img=bundle.images[1].pixels; // Anchor is not the first camera.
    MeasureGroup cached=bundle;cached.cloneImages();
    bundle.images[1].pixels.setTo(0);
    require(cached.img.data==cached.images[1].pixels.data && cached.img.at<uint8_t>(0,0)==9,
            "replay deep copy preserves a nonfirst anchor and isolates pixels");
    const double missing = std::numeric_limits<double>::infinity();
    require(cameraBundleReady(1,{1},{1},.0001,.5,false),"single-camera bundle");
    require(cameraBundleReady(1,{1,1.000025},{1,1.000025},.0001,.5,false),"25us synchronized views");
    require(!cameraBundleReady(1,{1,missing},{1,-1},.0001,.5,false),"wait for unresolved auxiliary image");
    require(cameraBundleReady(1,{1,1.1},{1,1.1},.0001,.5,false),"missing auxiliary closes on its own progress");
    require(cameraBundleReady(1,{1,missing},{1.6,-1},.0001,.5,false),"absent camera closes on sensor watermark");
    require(cameraBundleReady(1,{1,missing},{1,-1},.0001,.5,true),"EOF drains missing auxiliary");
    StatesGroup state;
    state.configureCameras(4);
    require(state.dimension() == 22, "four-camera dimension");
    state.exposure(3) = 1.2;
    state.cov(3, StatesGroup::exposureIndex(2)) = .00002;
    state.cov(StatesGroup::exposureIndex(2), 3) = .00002;
    StatesGroup prior = state;
    VectorXd delta = VectorXd::Zero(22);
    delta(3) = .05; delta(StatesGroup::exposureIndex(3)) = .1;
    state += delta;
    require(((state-prior)-delta).norm() < 1e-12, "state boxplus/boxminus preserves extra exposures");
    require(state.cov(3, StatesGroup::exposureIndex(2)) == .00002, "state copy preserves cross covariance");
    // A pose-only measurement must also update correlated auxiliary exposure.
    MatrixXd H = MatrixXd::Zero(1,22); H(0,3) = 1;
    MatrixXd gain = state.cov*H.transpose()*(H*state.cov*H.transpose()+MatrixXd::Identity(1,1)*.01).inverse();
    require(gain(StatesGroup::exposureIndex(2),0) != 0, "pose measurement updates exposure cross correlation");
    MatrixXd random = MatrixXd::Random(22,22);
    MatrixXd positive = random*random.transpose()+MatrixXd::Identity(22,22);
    require((inverseStateMatrix(positive)-positive.inverse()).norm()<1e-11,
            "block inverse equals full augmented covariance inverse");
    MatrixXd F = MatrixXd::Identity(22,22);
    F.block<3,3>(3,7) = M3D::Identity()*.01;
    MatrixXd Q = MatrixXd::Identity(22,22)*.0001;
    state.cov=positive;const MatrixXd expected=F*positive*F.transpose()+Q;
    propagateStateCovariance(state,F,Q);
    require((state.cov-expected).norm()<1e-11,"full covariance propagation preserves cross terms");

    unordered_map<VOXEL_LOCATION, VOXEL_POINTS *> map;
    for (size_t k = 0; k < 4; ++k)
    {
      VisualCamera view(map);
      view.camera_index = k; view.camera_count = 4; view.state = &prior;
      view.cam = new vk::PinholeCamera(64,48,1,30,31,32,24);
      view.grid_size = 12; view.grid_n_height = 4; view.patch_size = 4; view.patch_pyrimid_level = 1;
      view.raycast_en = false; view.normal_en = true; view.colmap_output_en = false;
      view.exposure_estimate_en = true; view.inverse_composition_en = false;
      view.setImuToLidarExtrinsic(V3D::Zero(), M3D::Identity());
      M3D r = Eigen::AngleAxisd(.05*k, V3D::UnitY()).toRotationMatrix();
      vector<double> rotation(9), translation={.02*double(k),.01*double(k),0};
      for (int i=0;i<3;++i) for(int j=0;j<3;++j) rotation[3*i+j]=r(i,j);
      view.setLidarToCameraExtrinsic(rotation,translation); view.initializeVIO();
      cv::Mat image(48,64,CV_8UC1);
      for(int y=0;y<48;++y) for(int x=0;x<64;++x) image.at<uint8_t>(y,x)=x+2*y;
      VisualPoint point(V3D(.1,.05,4),4);
      view.total_points=1;
      view.visual_submap->voxel_points={&point};
      view.visual_submap->search_levels={0};
      view.visual_submap->warp_patch={vector<float>(16,0)};
      view.visual_submap->inv_expo_list={1};
      view.visual_submap->errors={0};
      vector<double> weights={1};
      auto residual=[&](MatrixXd &jacobian) {
        jacobian=MatrixXd::Zero(16,7);VectorXd z=VectorXd::Zero(16);
        vector<float> errors(1); vector<int> counts(1);
        view.buildResiduals(image,0,weights,jacobian,z,errors,counts);
        require(counts[0]==16,"synthetic point contributes a complete patch");return z;
      };
      MatrixXd jacobian; residual(jacobian);
      for(int column=0;column<7;++column)
      {
        const double epsilon=.0005;
        VectorXd step=VectorXd::Zero(prior.dimension());
        const int index=column==6 ? StatesGroup::exposureIndex(k) : column;
        step(index)=epsilon;
        StatesGroup saved=prior;
        prior+=step;MatrixXd ignored;VectorXd plus=residual(ignored);
        prior=saved;prior+=(-step);VectorXd minus=residual(ignored);prior=saved;
        const VectorXd numerical=(plus-minus)/(2*epsilon);
        require((numerical-jacobian.col(column)).cwiseAbs().maxCoeff()<.04,
                "camera-specific pose/exposure photometric Jacobian finite difference");
      }
      view.visual_submap->reset();
      delete view.cam; view.cam=nullptr;
    }
    std::cout << "PASS: 22-state algebra/cross covariance and four extrinsic-specific photometric Jacobians\n";
  }
  catch (const std::exception &error)
  { std::cerr << error.what() << '\n'; return 1; }
  return 0;
}
