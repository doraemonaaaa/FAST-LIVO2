#ifndef FAST_LIVO_STATE_NUMERICS_H
#define FAST_LIVO_STATE_NUMERICS_H
#include "common_lib.h"

// Preserve the original fixed-size motion/primary-exposure arithmetic while
// appending camera states. These are block-matrix operations for all N, not
// separate single/multi-camera estimators. This also avoids changing Eigen's
// summation strategy for the existing block merely by making storage dynamic.
inline MatrixXd inverseStateMatrix(const MatrixXd &matrix)
{
  const MD(DIM_STATE,DIM_STATE) base = matrix.topLeftCorner<DIM_STATE,DIM_STATE>();
  const MD(DIM_STATE,DIM_STATE) inverse = base.inverse();
  MatrixXd result = MatrixXd::Zero(matrix.rows(), matrix.cols());
  result.topLeftCorner<DIM_STATE,DIM_STATE>() = inverse;
  const int extra = matrix.rows() - DIM_STATE;
  if (extra)
  {
    const MatrixXd B = matrix.topRightCorner(DIM_STATE,extra);
    const MatrixXd C = matrix.bottomLeftCorner(extra,DIM_STATE);
    const MatrixXd S = (matrix.bottomRightCorner(extra,extra)-C*inverse*B).inverse();
    result.topLeftCorner<DIM_STATE,DIM_STATE>() += inverse*B*S*C*inverse;
    result.topRightCorner(DIM_STATE,extra) = -inverse*B*S;
    result.bottomLeftCorner(extra,DIM_STATE) = -S*C*inverse;
    result.bottomRightCorner(extra,extra) = S;
  }
  return result;
}

inline void propagateStateCovariance(StatesGroup &state, const MatrixXd &F, const MatrixXd &Q)
{
  const MD(DIM_STATE,DIM_STATE) f = F.topLeftCorner<DIM_STATE,DIM_STATE>();
  const MD(DIM_STATE,DIM_STATE) p = state.cov.topLeftCorner<DIM_STATE,DIM_STATE>();
  const MD(DIM_STATE,DIM_STATE) q = Q.topLeftCorner<DIM_STATE,DIM_STATE>();
  const int extra = state.dimension()-DIM_STATE;
  state.cov.topLeftCorner<DIM_STATE,DIM_STATE>() = f*p*f.transpose()+q;
  if (extra)
  {
    state.cov.topRightCorner(DIM_STATE,extra) = (f*state.cov.topRightCorner(DIM_STATE,extra)).eval();
    state.cov.bottomLeftCorner(extra,DIM_STATE) = (state.cov.bottomLeftCorner(extra,DIM_STATE)*f.transpose()).eval();
    state.cov.bottomRightCorner(extra,extra) += Q.bottomRightCorner(extra,extra);
  }
}

template<int Columns>
inline MatrixXd gainInformation(const MatrixXd &K, const Matrix<double,Columns,Columns> &information)
{
  MatrixXd result(K.rows(), Columns);
  result.topRows<DIM_STATE>() = K.block<DIM_STATE,Columns>(0,0)*information;
  const int extra = K.rows()-DIM_STATE;
  if (extra) result.bottomRows(extra) = K.block(DIM_STATE,0,extra,Columns)*information;
  return result;
}
#endif
