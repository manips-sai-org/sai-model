#ifndef RBDL_SPATIALALGEBRAOPERATORS_ED_H
#define RBDL_SPATIALALGEBRAOPERATORS_ED_H

#include <rbdl/rbdl.h>
#include <rbdl/SpatialAlgebraOperators.h>

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace Math {
// -----------------------------------------------------------------------------
namespace ED {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;

typedef Eigen::Matrix<double, 6, Eigen::Dynamic> SpatialDirection;
typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Direction3d;

typedef Eigen::Ref<Eigen::Vector3d> Vector3d_ref;
typedef Eigen::Ref<const Eigen::Vector3d> Vector3d_cref;

typedef Eigen::Ref<Direction3d> Direction3d_ref;
typedef Eigen::Ref<const Direction3d> Direction3d_cref;

/*
inline void X_apply_v_plus_w (
  SpatialVector &res,
  SpatialDirection &res_dir,
  const SpatialTransform &X, const std::vector<SpatialTransform> &X_dir,
  const SpatialVector &v, const SpatialDirection &v_dir,
  const SpatialVector &w
) {
  Vector3d v_rxw (
      v[3] - X.r[1]*v[2] + X.r[2]*v[1],
      v[4] - X.r[2]*v[0] + X.r[0]*v[2],
      v[5] - X.r[0]*v[1] + X.r[1]*v[0]
  );

  const unsigned int ndirs = res_dir.cols();
  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d v_rxw_dir (
        v_dir(3, idir) - X.r[1]*v_dir(2, idir) + X.r[2]*v_dir(1, idir),
        v_dir(4, idir) - X.r[2]*v_dir(0, idir) + X.r[0]*v_dir(2, idir),
        v_dir(5, idir) - X.r[0]*v_dir(1, idir) + X.r[1]*v_dir(0, idir)
    );

      res_dir(0, idir) = X_dir[idir].E(0,0) * v[0]     + X_dir[idir].E(0,1) * v[1]     + X_dir[idir].E(0,2) * v[2]
                       + X.E(0,0)           * v_dir(0, idir) + X.E(0,1)     * v_dir(1, idir) + X.E(0,2)     * v_dir(2, idir);

      res_dir(1, idir) = X_dir[idir].E(1,0) * v[0]     + X_dir[idir].E(1,1) * v[1]     + X_dir[idir].E(1,2) * v[2]
             + X.E(1,0)     * v_dir(0, idir) + X.E(1,1)     * v_dir(1, idir) + X.E(1,2)     * v_dir(2, idir);

      res_dir(2, idir) = X_dir[idir].E(2,0) * v[0]     + X_dir[idir].E(2,1) * v[1]     + X_dir[idir].E(2,2) * v[2]
             + X.E(2,0)     * v_dir(0, idir) + X.E(2,1)     * v_dir(1, idir) + X.E(2,2)     * v_dir(2, idir);

      res_dir(3, idir) = X_dir[idir].E(0,0) * v_rxw[0]     + X_dir[idir].E(0,1) * v_rxw[1]     + X_dir[idir].E(0,2) * v_rxw[2]
             + X.E(0,0)     * v_rxw_dir[0] + X.E(0,1)     * v_rxw_dir[1] + X.E(0,2)     * v_rxw_dir[2];

      res_dir(4, idir) = X_dir[idir].E(1,0) * v_rxw[0]     + X_dir[idir].E(1,1) * v_rxw[1]     + X_dir[idir].E(1,2) * v_rxw[2]
             + X.E(1,0)     * v_rxw_dir[0] + X.E(1,1)     * v_rxw_dir[1] + X.E(1,2)     * v_rxw_dir[2];

      res_dir(5, idir) = X_dir[idir].E(2,0) * v_rxw[0]     + X_dir[idir].E(2,1) * v_rxw[1]     + X_dir[idir].E(2,2) * v_rxw[2]
             + X.E(2,0)     * v_rxw_dir[0] + X.E(2,1)     * v_rxw_dir[1] + X.E(2,2)     * v_rxw_dir[2];
    }

    // nominal evaluation
    res[0] = X.E(0,0) * v[0] + X.E(0,1) * v[1] + X.E(0,2) * v[2];
    res[1] = X.E(1,0) * v[0] + X.E(1,1) * v[1] + X.E(1,2) * v[2];
    res[2] = X.E(2,0) * v[0] + X.E(2,1) * v[1] + X.E(2,2) * v[2];
    res[3] = X.E(0,0) * v_rxw[0] + X.E(0,1) * v_rxw[1] + X.E(0,2) * v_rxw[2];
    res[4] = X.E(1,0) * v_rxw[0] + X.E(1,1) * v_rxw[1] + X.E(1,2) * v_rxw[2];
    res[5] = X.E(2,0) * v_rxw[0] + X.E(2,1) * v_rxw[1] + X.E(2,2) * v_rxw[2];
  res = SpatialVector(
      X.E(0,0) * v[0] + X.E(0,1) * v[1] + X.E(0,2) * v[2] + u[0],
      X.E(1,0) * v[0] + X.E(1,1) * v[1] + X.E(1,2) * v[2] + u[1],
      X.E(2,0) * v[0] + X.E(2,1) * v[1] + X.E(2,2) * v[2] + u[2],
      X.E(0,0) * v_rxw[0] + X.E(0,1) * v_rxw[1] + X.E(0,2) * v_rxw[2] + u[3],
      X.E(1,0) * v_rxw[0] + X.E(1,1) * v_rxw[1] + X.E(1,2) * v_rxw[2] + u[4],
      X.E(2,0) * v_rxw[0] + X.E(2,1) * v_rxw[1] + X.E(2,2) * v_rxw[2] + u[5]
      );
}
*/

/** Same as res = X_d * v + X * v_d.
 *
 * \returns (E * w, - E * rxw + E * v)
 */
inline void X_apply_v_add_u (
  SpatialVector &res,
  SpatialDirection &res_dir,
  const SpatialTransform &X, const std::vector<SpatialTransform> &X_dir,
  const SpatialVector &v, const SpatialDirection &v_dir,
  const SpatialVector &u, const SpatialDirection &u_dir,
  const unsigned int &ndirs
) {

  Vector3d v_rxw (
      v[3] - X.r[1] * v[2] + X.r[2] * v[1],
      v[4] - X.r[2] * v[0] + X.r[0] * v[2],
      v[5] - X.r[0] * v[1] + X.r[1] * v[0]
      );

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d w_dir = v_dir.col(idir).segment<3>(0);
    res_dir.col(idir).segment<3>(0) =
        X.E * w_dir
        + X_dir[idir].E * v.segment<3>(0);
    res_dir.col(idir).segment<3>(3) =
        X.E * (v_dir.col(idir).segment<3>(3)
                - X.r.cross(w_dir)
                - X_dir[idir].r.cross(v.segment<3>(0)))
        + X_dir[idir].E * v_rxw;
  }
  res_dir.leftCols(ndirs) += u_dir.leftCols(ndirs);

  res = SpatialVector(
      X.E(0,0) * v[0] + X.E(0,1) * v[1] + X.E(0,2) * v[2] + u[0],
      X.E(1,0) * v[0] + X.E(1,1) * v[1] + X.E(1,2) * v[2] + u[1],
      X.E(2,0) * v[0] + X.E(2,1) * v[1] + X.E(2,2) * v[2] + u[2],
      X.E(0,0) * v_rxw[0] + X.E(0,1) * v_rxw[1] + X.E(0,2) * v_rxw[2] + u[3],
      X.E(1,0) * v_rxw[0] + X.E(1,1) * v_rxw[1] + X.E(1,2) * v_rxw[2] + u[4],
      X.E(2,0) * v_rxw[0] + X.E(2,1) * v_rxw[1] + X.E(2,2) * v_rxw[2] + u[5]
      );
}

/** Same as res = X_d * v + X * v_d.
 *
 * \returns (E * w, - E * rxw + E * v)
 */
inline void X_apply_v (
  SpatialVector &res,
  SpatialDirection &res_dir,
  const SpatialTransform &X, const std::vector<SpatialTransform> &X_dir,
  const SpatialVector &v, const SpatialDirection &v_dir,
  const unsigned int &ndirs
) {

  Vector3d v_rxw (
      v[3] - X.r[1] * v[2] + X.r[2] * v[1],
      v[4] - X.r[2] * v[0] + X.r[0] * v[2],
      v[5] - X.r[0] * v[1] + X.r[1] * v[0]
      );

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d w_dir = v_dir.col(idir).segment<3>(0);
    res_dir.col(idir).segment<3>(0) =
        X.E * w_dir
        + X_dir[idir].E * v.segment<3>(0);
    res_dir.col(idir).segment<3>(3) =
        X.E * (v_dir.col(idir).segment<3>(3)
                - X.r.cross(w_dir)
                - X_dir[idir].r.cross(v.segment<3>(0)))
        + X_dir[idir].E * v_rxw;
  }

  res = SpatialVector(
      X.E(0,0) * v[0] + X.E(0,1) * v[1] + X.E(0,2) * v[2],
      X.E(1,0) * v[0] + X.E(1,1) * v[1] + X.E(1,2) * v[2],
      X.E(2,0) * v[0] + X.E(2,1) * v[1] + X.E(2,2) * v[2],
      X.E(0,0) * v_rxw[0] + X.E(0,1) * v_rxw[1] + X.E(0,2) * v_rxw[2],
      X.E(1,0) * v_rxw[0] + X.E(1,1) * v_rxw[1] + X.E(1,2) * v_rxw[2],
      X.E(2,0) * v_rxw[0] + X.E(2,1) * v_rxw[1] + X.E(2,2) * v_rxw[2]
      );
}

/** Same as f += X^T * f.
 */
inline void inplace_X_applyTranspose_f(
  SpatialVector &res,
  SpatialDirection &res_dir,
  const SpatialTransform &X, const std::vector<SpatialTransform> &X_dir,
  const SpatialVector &f_sp, const SpatialDirection &f_sp_dir,
  const unsigned int &ndirs
) {
  Vector3d E_T_f (
    X.E(0,0) * f_sp[3] + X.E(1,0) * f_sp[4] + X.E(2,0) * f_sp[5],
    X.E(0,1) * f_sp[3] + X.E(1,1) * f_sp[4] + X.E(2,1) * f_sp[5],
    X.E(0,2) * f_sp[3] + X.E(1,2) * f_sp[4] + X.E(2,2) * f_sp[5]
  );

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d E_T_f_dir (
      X_dir[idir].E(0,0) * f_sp[3] + X_dir[idir].E(1,0) * f_sp[4] + X_dir[idir].E(2,0) * f_sp[5] + X.E(0,0) * f_sp_dir(3, idir) + X.E(1,0) * f_sp_dir(4, idir) + X.E(2,0) * f_sp_dir(5, idir),
      X_dir[idir].E(0,1) * f_sp[3] + X_dir[idir].E(1,1) * f_sp[4] + X_dir[idir].E(2,1) * f_sp[5] + X.E(0,1) * f_sp_dir(3, idir) + X.E(1,1) * f_sp_dir(4, idir) + X.E(2,1) * f_sp_dir(5, idir),
      X_dir[idir].E(0,2) * f_sp[3] + X_dir[idir].E(1,2) * f_sp[4] + X_dir[idir].E(2,2) * f_sp[5] + X.E(0,2) * f_sp_dir(3, idir) + X.E(1,2) * f_sp_dir(4, idir) + X.E(2,2) * f_sp_dir(5, idir)
    );

    res_dir(0, idir) += X_dir[idir].E(0,0) * f_sp[0]     + X_dir[idir].E(1,0) * f_sp[1]     + X_dir[idir].E(2,0) * f_sp[2]     - X_dir[idir].r[2] * E_T_f[1]     + X_dir[idir].r[1] * E_T_f[2]
           + X.E(0,0)     * f_sp_dir(0, idir) + X.E(1,0)     * f_sp_dir(1, idir) + X.E(2,0)     * f_sp_dir(2, idir) - X.r[2]     * E_T_f_dir[1] + X.r[1]     * E_T_f_dir[2];
    res_dir(1, idir) += X_dir[idir].E(0,1) * f_sp[0]     + X_dir[idir].E(1,1) * f_sp[1]     + X_dir[idir].E(2,1) * f_sp[2]     + X_dir[idir].r[2] * E_T_f[0]     - X_dir[idir].r[0] * E_T_f[2]
           + X.E(0,1)     * f_sp_dir(0, idir) + X.E(1,1)     * f_sp_dir(1, idir) + X.E(2,1)     * f_sp_dir(2, idir) + X.r[2]     * E_T_f_dir[0] - X.r[0]     * E_T_f_dir[2];
    res_dir(2, idir) += X_dir[idir].E(0,2) * f_sp[0]     + X_dir[idir].E(1,2) * f_sp[1]     + X_dir[idir].E(2,2) * f_sp[2]     - X_dir[idir].r[1] * E_T_f[0]     + X_dir[idir].r[0] * E_T_f[1]
           + X.E(0,2)     * f_sp_dir(0, idir) + X.E(1,2)     * f_sp_dir(1, idir) + X.E(2,2)     * f_sp_dir(2, idir) - X.r[1]     * E_T_f_dir[0] + X.r[0]     * E_T_f_dir[1];
    res_dir(3, idir) += E_T_f_dir[0];
    res_dir(4, idir) += E_T_f_dir[1];
    res_dir(5, idir) += E_T_f_dir[2];
  }

  res[0] += X.E(0,0) * f_sp[0] + X.E(1,0) * f_sp[1] + X.E(2,0) * f_sp[2] - X.r[2] * E_T_f[1] + X.r[1] * E_T_f[2];
  res[1] += X.E(0,1) * f_sp[0] + X.E(1,1) * f_sp[1] + X.E(2,1) * f_sp[2] + X.r[2] * E_T_f[0] - X.r[0] * E_T_f[2];
  res[2] += X.E(0,2) * f_sp[0] + X.E(1,2) * f_sp[1] + X.E(2,2) * f_sp[2] - X.r[1] * E_T_f[0] + X.r[0] * E_T_f[1];
  res[3] += E_T_f [0];
  res[4] += E_T_f [1];
  res[5] += E_T_f [2];

  return;
}

/** Same as X^T * f.
 */
inline void X_applyT_f(
  SpatialVector &res,
  SpatialDirection &res_dir,
  const SpatialTransform &X, const std::vector<SpatialTransform> &X_dir,
  const SpatialVector &f_sp, const SpatialDirection &f_sp_dir,
  const unsigned int &ndirs
) {
  Vector3d E_T_f (
    X.E(0,0) * f_sp[3] + X.E(1,0) * f_sp[4] + X.E(2,0) * f_sp[5],
    X.E(0,1) * f_sp[3] + X.E(1,1) * f_sp[4] + X.E(2,1) * f_sp[5],
    X.E(0,2) * f_sp[3] + X.E(1,2) * f_sp[4] + X.E(2,2) * f_sp[5]
  );

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d E_T_f_dir (
      X_dir[idir].E(0,0) * f_sp[3] + X_dir[idir].E(1,0) * f_sp[4] + X_dir[idir].E(2,0) * f_sp[5] + X.E(0,0) * f_sp_dir(3, idir) + X.E(1,0) * f_sp_dir(4, idir) + X.E(2,0) * f_sp_dir(5, idir),
      X_dir[idir].E(0,1) * f_sp[3] + X_dir[idir].E(1,1) * f_sp[4] + X_dir[idir].E(2,1) * f_sp[5] + X.E(0,1) * f_sp_dir(3, idir) + X.E(1,1) * f_sp_dir(4, idir) + X.E(2,1) * f_sp_dir(5, idir),
      X_dir[idir].E(0,2) * f_sp[3] + X_dir[idir].E(1,2) * f_sp[4] + X_dir[idir].E(2,2) * f_sp[5] + X.E(0,2) * f_sp_dir(3, idir) + X.E(1,2) * f_sp_dir(4, idir) + X.E(2,2) * f_sp_dir(5, idir)
    );

    res_dir(0, idir) = X_dir[idir].E(0,0) * f_sp[0]     + X_dir[idir].E(1,0) * f_sp[1]     + X_dir[idir].E(2,0) * f_sp[2]     - X_dir[idir].r[2] * E_T_f[1]     + X_dir[idir].r[1] * E_T_f[2]
           + X.E(0,0)     * f_sp_dir(0, idir) + X.E(1,0)     * f_sp_dir(1, idir) + X.E(2,0)     * f_sp_dir(2, idir) - X.r[2]     * E_T_f_dir[1] + X.r[1]     * E_T_f_dir[2];
    res_dir(1, idir) = X_dir[idir].E(0,1) * f_sp[0]     + X_dir[idir].E(1,1) * f_sp[1]     + X_dir[idir].E(2,1) * f_sp[2]     + X_dir[idir].r[2] * E_T_f[0]     - X_dir[idir].r[0] * E_T_f[2]
           + X.E(0,1)     * f_sp_dir(0, idir) + X.E(1,1)     * f_sp_dir(1, idir) + X.E(2,1)     * f_sp_dir(2, idir) + X.r[2]     * E_T_f_dir[0] - X.r[0]     * E_T_f_dir[2];
    res_dir(2, idir) = X_dir[idir].E(0,2) * f_sp[0]     + X_dir[idir].E(1,2) * f_sp[1]     + X_dir[idir].E(2,2) * f_sp[2]     - X_dir[idir].r[1] * E_T_f[0]     + X_dir[idir].r[0] * E_T_f[1]
           + X.E(0,2)     * f_sp_dir(0, idir) + X.E(1,2)     * f_sp_dir(1, idir) + X.E(2,2)     * f_sp_dir(2, idir) - X.r[1]     * E_T_f_dir[0] + X.r[0]     * E_T_f_dir[1];
    res_dir(3, idir) = E_T_f_dir[0];
    res_dir(4, idir) = E_T_f_dir[1];
    res_dir(5, idir) = E_T_f_dir[2];
  }

  res[0] = X.E(0,0) * f_sp[0] + X.E(1,0) * f_sp[1] + X.E(2,0) * f_sp[2] - X.r[2] * E_T_f[1] + X.r[1] * E_T_f[2];
  res[1] = X.E(0,1) * f_sp[0] + X.E(1,1) * f_sp[1] + X.E(2,1) * f_sp[2] + X.r[2] * E_T_f[0] - X.r[0] * E_T_f[2];
  res[2] = X.E(0,2) * f_sp[0] + X.E(1,2) * f_sp[1] + X.E(2,2) * f_sp[2] - X.r[1] * E_T_f[0] + X.r[0] * E_T_f[1];
  res[3] = E_T_f [0];
  res[4] = E_T_f [1];
  res[5] = E_T_f [2];

  return;
}


/** Same as X^* I X^{-1}
*/
// void apply (
//   SpatialRigidBodyInertia & res,
//   const SpatialTransform &X,
//   const SpatialTransform &X_dir,
//   const SpatialRigidBodyInertia &rbi,
//   const SpatialRigidBodyInertia &rbi_dir
// ) {
//   std::cerr << __FILE__ << " " << __LINE__ << ":"
//        << " not yet implemented." << std::endl;
//   abort();
//   // return SpatialRigidBodyInertia (
//   //     rbi.m,
//   //     E * (rbi.h - rbi.m * r),
//   //     E *
//   //     (
//   //      Matrix3d (
//   //        rbi.Ixx, rbi.Iyx, rbi.Izx,
//   //        rbi.Iyx, rbi.Iyy, rbi.Izy,
//   //        rbi.Izx, rbi.Izy, rbi.Izz
//   //        )
//   //      + VectorCrossMatrix (r) * VectorCrossMatrix (rbi.h)
//   //      + (VectorCrossMatrix(rbi.h - rbi.m * r) * VectorCrossMatrix (r))
//   //     )
//   //     * E.transpose()
//   //     );
// }


inline void Xrotx (
  SpatialTransform &res,
  std::vector<SpatialTransform> &res_dirs,
  const double &rot,
  const VectorNd &rot_dirs,
  const SpatialTransform &X,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (rot);
  c = cos (rot);

  res.E = Matrix3d (
    1., 0., 0.,
    0., c, s,
    0., -s, c
  );
  res.r = X.r;

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*rot_dirs[idir];
    c_t = c*rot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
        0.,   0.,   0.,
        0., -s_t,  c_t,
        0., -c_t, -s_t
    ) * X.E;
    res_dirs[idir].r.setZero();
  }

  return;
}

inline void Xrotx (
  SpatialTransform &res,
  std::vector<SpatialTransform> &res_dirs,
  const double &rot,
  const VectorNd &rot_dirs,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (rot);
  c = cos (rot);

  res.E = Matrix3d (
    1., 0., 0.,
    0., c, s,
    0., -s, c
  );
  res.r.setZero();

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*rot_dirs[idir];
    c_t = c*rot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
        0.,   0.,   0.,
        0., -s_t,  c_t,
        0., -c_t, -s_t
    );
    res_dirs[idir].r.setZero();
  }

  return;
}

inline void Xroty (
  SpatialTransform &res,
  std::vector<SpatialTransform> &res_dirs,
  const double &rot,
  const VectorNd &rot_dirs,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (rot);
  c = cos (rot);

  res.E = Matrix3d (
    c, 0., -s,
    0., 1., 0.,
    s, 0., c
  );
  res.r.setZero();

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*rot_dirs[idir];
    c_t = c*rot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
     -s_t, 0., -c_t,
       0., 0.,   0.,
      c_t, 0., -s_t
    );
    res_dirs[idir].r.setZero();
  }

  return;
}

inline void Xroty (
  SpatialTransform &res,
  std::vector<SpatialTransform> &res_dirs,
  const double &rot,
  const VectorNd &rot_dirs,
  const SpatialTransform &X,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (rot);
  c = cos (rot);

  res.E = Matrix3d (
    c, 0., -s,
    0., 1., 0.,
    s, 0., c
  );
  res.r = X.r; // + X.E.transpose() * Vector3d (0., 0., 0.);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*rot_dirs[idir];
    c_t = c*rot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
     -s_t, 0., -c_t,
       0., 0.,   0.,
      c_t, 0., -s_t
    ) * X.E;
    res_dirs[idir].r.setZero();
  }

  return;
}


inline void Xrotz (
  SpatialTransform &res,
  std::vector<SpatialTransform> &res_dirs,
  const double &rot,
  const VectorNd &rot_dirs,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (rot);
  c = cos (rot);

  res.E = Matrix3d (
    c, s, 0.,
    -s, c, 0.,
    0., 0., 1.
  );
  res.r.setZero();

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*rot_dirs[idir];
    c_t = c*rot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
      -s_t,  c_t, 0.,
      -c_t, -s_t, 0.,
        0.,   0., 0.
    );
    res_dirs[idir].r.setZero();
  }

  return;
}

inline void Xrotz (
  SpatialTransform &res,
  std::vector<SpatialTransform> &res_dirs,
  const double &rot,
  const VectorNd &rot_dirs,
  const SpatialTransform &X,
  const unsigned int &ndirs
) {
  double s, c, s_t, c_t;
  // NOTE sin and cos are slow to compute?
  s = sin (rot);
  c = cos (rot);

  res.E = Matrix3d (
    c, s, 0.,
    -s, c, 0.,
    0., 0., 1.
  );
  res.r = X.r; // + X.E.transpose() * Vector3d (0., 0., 0.);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    s_t = s*rot_dirs[idir];
    c_t = c*rot_dirs[idir];
    res_dirs[idir].E = Matrix3d (
      -s_t,  c_t, 0.,
      -c_t, -s_t, 0.,
        0.,   0., 0.
    ) * X.E;
    res_dirs[idir].r.setZero();
  }

  return;
}

inline void crossm (
  SpatialVector &res, SpatialDirection &res_dir,
  const SpatialVector &v1, const SpatialDirection &v1_dir,
  const SpatialVector &v2, const SpatialDirection &v2_dir,
  const unsigned int &ndirs
) {
  // NOTE we require row dimension to be 3 at compile time, we use Eigen::Ref
  //      to achieve this
  // NOTE we use identity a x b = - b x a and use colwise partial reduction
  //      to efficiently perform vectorized cross product evaluations
  const Vector3d_cref v1_head = v1.head(3);
  const Vector3d_cref v2_head = v2.head(3);
  const Vector3d_cref v1_tail = v1.tail(3);
  const Vector3d_cref v2_tail = v2.tail(3);

  const Direction3d_cref res_dir_head = res_dir.leftCols(ndirs).topRows(3);
  const Direction3d_cref v1_dir_head  = v1_dir.leftCols(ndirs).topRows(3);
  const Direction3d_cref v2_dir_head  = v2_dir.leftCols(ndirs).topRows(3);

  const Direction3d_cref res_dir_tail = res_dir.leftCols(ndirs).bottomRows(3);
  const Direction3d_cref v1_dir_tail  = v1_dir.leftCols(ndirs).bottomRows(3);
  const Direction3d_cref v2_dir_tail  = v2_dir.leftCols(ndirs).bottomRows(3);

  res_dir.leftCols(ndirs).topRows(3)
    = v1_dir_head.colwise().cross(v2_head)
    - v2_dir_head.colwise().cross(v1_head);

  res_dir.leftCols(ndirs).bottomRows(3)
    = v1_dir_tail.colwise().cross(v2_head)
    - v2_dir_head.colwise().cross(v1_tail)
    + v1_dir_head.colwise().cross(v2_tail)
    - v2_dir_tail.colwise().cross(v1_head);

  res.head(3) = v1_head.cross(v2_head);
  res.tail(3) = v1_tail.cross(v2_head) + v1_head.cross(v2_tail);

  return;
}

inline void crossf (
  SpatialVector &res, SpatialDirection &res_dir,
  const SpatialVector &v1, const SpatialDirection &v1_dir,
  const SpatialVector &v2, const SpatialDirection &v2_dir,
  const unsigned int &ndirs
) {
  // NOTE we require row dimension to be 3 at compile time, we use Eigen::Ref
  //      to achieve this
  // NOTE we use identity a x b = - b x a and use colwise partial reduction
  //      to efficiently perform vectorized cross product evaluations
  const Vector3d_cref v1_head = v1.head(3);
  const Vector3d_cref v2_head = v2.head(3);
  const Vector3d_cref v1_tail = v1.tail(3);
  const Vector3d_cref v2_tail = v2.tail(3);

  const Direction3d_cref v1_dir_head = v1_dir.leftCols(ndirs).topRows(3);
  const Direction3d_cref v2_dir_head = v2_dir.leftCols(ndirs).topRows(3);

  const Direction3d_cref v1_dir_tail = v1_dir.leftCols(ndirs).bottomRows(3);
  const Direction3d_cref v2_dir_tail = v2_dir.leftCols(ndirs).bottomRows(3);

  res_dir.leftCols(ndirs).topRows(3)
    = v1_dir_head.colwise().cross(v2_head)
    - v2_dir_head.colwise().cross(v1_head)
    + v1_dir_tail.colwise().cross(v2_tail)
    - v2_dir_tail.colwise().cross(v1_tail);

  res_dir.leftCols(ndirs).bottomRows(3)
    = v1_dir_head.colwise().cross(v2_tail)
    - v2_dir_tail.colwise().cross(v1_head);

  res.segment<3>(0) = v1_head.cross(v2_head) + v1_tail.cross(v2_tail);
  res.segment<3>(3) = v1_head.cross(v2_tail);

  return;
}


inline void SRBI_apply_v (
  SpatialVector &res, SpatialDirection &res_dir,
  const SpatialRigidBodyInertia &I,
  const SpatialVector &v, const SpatialDirection &v_dir,
  const unsigned int &ndirs
) {
  const Vector3d_cref v_head = v.head(3);
  const Vector3d_cref v_tail = v.tail(3);
  const Direction3d_cref v_dir_head = v_dir.leftCols(ndirs).topRows(3);
  const Direction3d_cref v_dir_tail = v_dir.leftCols(ndirs).bottomRows(3);

  Matrix3d I_m = Matrix3d (
      I.Ixx, I.Iyx, I.Izx,
      I.Iyx, I.Iyy, I.Izy,
      I.Izx, I.Izy, I.Izz
  );

  // derivative evaluation
  res_dir.leftCols(ndirs).topRows(3)
    = I_m * v_dir_head - v_dir_tail.colwise().cross(I.h);

  res_dir.leftCols(ndirs).bottomRows(3)
    = I.m * v_dir_tail - v_dir_head.colwise().cross(I.h);

  // nominal evaluation
  res.head(3) = I_m * v_head + I.h.cross(v_tail);
  res.tail(3) = I.m * v_tail - I.h.cross(v_head);

  return;
}

inline void add_SRBI_apply_v  (
  SpatialVector &res, SpatialDirection &res_dir,
  const SpatialRigidBodyInertia &I,
  const SpatialVector &v, const SpatialDirection &v_dir,
  const unsigned int &ndirs
) {
  const Vector3d_cref v_head = v.head(3);
  const Vector3d_cref v_tail = v.tail(3);
  const Direction3d_cref v_dir_head = v_dir.leftCols(ndirs).topRows(3);
  const Direction3d_cref v_dir_tail = v_dir.leftCols(ndirs).bottomRows(3);

  Matrix3d I_m = Matrix3d (
      I.Ixx, I.Iyx, I.Izx,
      I.Iyx, I.Iyy, I.Izy,
      I.Izx, I.Izy, I.Izz
  );

  // derivative evaluation
  res_dir.leftCols(ndirs).topRows(3)
    += I_m * v_dir_head - v_dir_tail.colwise().cross(I.h);

  res_dir.leftCols(ndirs).bottomRows(3)
    += I.m * v_dir_tail - v_dir_head.colwise().cross(I.h);

  // nominal evaluation
  res.head(3) += I_m * v_head + I.h.cross(v_tail);
  res.tail(3) += I.m * v_tail - I.h.cross(v_head);

  return;
}

// -----------------------------------------------------------------------------
} // namespace ED
// -----------------------------------------------------------------------------
} // namespace Math
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif
