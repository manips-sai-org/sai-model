#ifndef RBDL_SPATIALALGEBRAOPERATORS_AD_H
#define RBDL_SPATIALALGEBRAOPERATORS_AD_H

#include <rbdl/rbdl.h>
#include <rbdl/SpatialAlgebraOperators.h>

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;

/**
 * @brief mulSTST for ST A and ST B compute A * B
 *        save against &lhs == &res, &rhs == &res,
 *                     &res_dirs == &lhs_dirs, &res_dirs == &rhs_dirs
 * @param lhs A
 * @param lhs_dirs dA
 * @param rhs B
 * @param rhs_dirs
 * @param res A * B
 * @param res_dirs d(A * B)
 */
inline void mulSTST(
    unsigned ndirs,
    SpatialTransform const & lhs,
    std::vector<SpatialTransform> const & lhs_dirs,
    SpatialTransform const & rhs,
    std::vector<SpatialTransform> const & rhs_dirs,
    SpatialTransform & res,
    std::vector<SpatialTransform> & res_dirs) {
  assert(ndirs <= lhs_dirs.size());
  assert(ndirs <= rhs_dirs.size());
  assert(ndirs <= res_dirs.size());

  // derivative code
  for (unsigned idir = 0; idir < ndirs; idir++) {
    res_dirs[idir].r =
        rhs_dirs[idir].r
        + rhs.E.transpose() * lhs_dirs[idir].r
        + rhs_dirs[idir].E.transpose() * lhs.r;
    res_dirs[idir].E = lhs.E * rhs_dirs[idir].E + lhs_dirs[idir].E * rhs.E;
  }

  // nominal code
  res.r = rhs.r + rhs.E.transpose() * lhs.r;
  res.E = lhs.E * rhs.E;
}

/**
 * @brief mulSTST for ST A and ST B compute A * B if B is constant (dB = 0),
 *        save against &lhs == &res, &rhs == &res, &res_dirs == &lhs_dirs
 * @param lhs A
 * @param lhs_dirs dA
 * @param rhs B
 * @param res A * B
 * @param res_dirs d(A * B)
 */
inline void mulSTST(
    unsigned ndirs,
    SpatialTransform const & lhs,
    std::vector<SpatialTransform> const & lhs_dirs,
    SpatialTransform const & rhs,
    SpatialTransform & res,
    std::vector<SpatialTransform> & res_dirs) {
  assert(ndirs <= lhs_dirs.size());
  assert(ndirs <= res_dirs.size());

  // derivative code
  for (unsigned idir = 0; idir < ndirs; idir++) {
    res_dirs[idir].r = rhs.E.transpose() * lhs_dirs[idir].r;
    res_dirs[idir].E = lhs_dirs[idir].E * rhs.E;
  }

  // nominal code
  res.E = lhs.E * rhs.E;
  res.r = rhs.r + rhs.E.transpose() * lhs.r;
}

/**
 * @brief add_sqrFormSTSM_noalias ST A, SM B and SM C compute C += A^T * B * A
 * @param st A
 * @param st_dirs dA
 * @param sm B
 * @param sm_dirs dB
 * @param res C += A^T * B * A
 * @param res_dirs dC += d(A^T * B * A)
 */
inline void addSqrFormSTSM_noalias (
    unsigned ndirs,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialMatrix const & sm,
    std::vector<SpatialMatrix> const & sm_dirs,
    SpatialMatrix & res,
    std::vector<SpatialMatrix> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sm_dirs.size());
  assert(ndirs <= res_dirs.size());

  SpatialMatrix A   = st.toMatrix();
  SpatialMatrix A_T = A.transpose();

  for (unsigned idir = 0; idir < ndirs; idir++) {
    SpatialMatrix A_dir;
    A_dir.block<3,3>(0,0) = st_dirs[idir].E;
    A_dir.block<3,3>(3,0) =
        - st_dirs[idir].E * VectorCrossMatrix(st.r)
        - st.E * VectorCrossMatrix(st_dirs[idir].r);
    A_dir.block<3,3>(0,3).setZero();
    A_dir.block<3,3>(3,3) = st_dirs[idir].E;

    res_dirs[idir].noalias() +=
        A_dir.transpose() * sm * A
        + A_T * sm_dirs[idir] * A
        +  A_T * sm * A_dir;
  }
  res.noalias() += A_T * sm * A;
}

inline void applySTSV (
    unsigned ndirs,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    SpatialVector & res,
    std::vector<SpatialVector> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= res_dirs.size());

  Vector3d v_rxw (
      sv[3] - st.r[1] * sv[2] + st.r[2] * sv[1],
      sv[4] - st.r[2] * sv[0] + st.r[0] * sv[2],
      sv[5] - st.r[0] * sv[1] + st.r[1] * sv[0]
      );

  for (unsigned idir = 0; idir < ndirs; idir++) {
    res_dirs[idir].segment<3>(0) =
        st_dirs[idir].E * sv.segment<3>(0);
    res_dirs[idir].segment<3>(3) =
        st.E * (- st_dirs[idir].r.cross(sv.segment<3>(0)))
        + st_dirs[idir].E * v_rxw;
  }

  res = SpatialVector(
      st.E(0,0) * sv[0] + st.E(0,1) * sv[1] + st.E(0,2) * sv[2],
      st.E(1,0) * sv[0] + st.E(1,1) * sv[1] + st.E(1,2) * sv[2],
      st.E(2,0) * sv[0] + st.E(2,1) * sv[1] + st.E(2,2) * sv[2],
      st.E(0,0) * v_rxw[0] + st.E(0,1) * v_rxw[1] + st.E(0,2) * v_rxw[2],
      st.E(1,0) * v_rxw[0] + st.E(1,1) * v_rxw[1] + st.E(1,2) * v_rxw[2],
      st.E(2,0) * v_rxw[0] + st.E(2,1) * v_rxw[1] + st.E(2,2) * v_rxw[2]
      );
}

inline void applySTSV(
    unsigned ndirs,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    std::vector<SpatialVector> const & sv_dirs,
    SpatialVector & res,
    std::vector<SpatialVector> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.size());
  assert(ndirs <= res_dirs.size());

  Vector3d v_rxw (
      sv[3] - st.r[1] * sv[2] + st.r[2] * sv[1],
      sv[4] - st.r[2] * sv[0] + st.r[0] * sv[2],
      sv[5] - st.r[0] * sv[1] + st.r[1] * sv[0]
      );

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d w_dir = sv_dirs[idir].segment<3>(0);
    res_dirs[idir].segment<3>(0) =
        st.E * w_dir
        + st_dirs[idir].E * sv.segment<3>(0);
    res_dirs[idir].segment<3>(3) =
        st.E * (sv_dirs[idir].segment<3>(3)
                - st.r.cross(w_dir)
                - st_dirs[idir].r.cross(sv.segment<3>(0)))
        + st_dirs[idir].E * v_rxw;
  }

  res = SpatialVector(
      st.E(0,0) * sv[0] + st.E(0,1) * sv[1] + st.E(0,2) * sv[2],
      st.E(1,0) * sv[0] + st.E(1,1) * sv[1] + st.E(1,2) * sv[2],
      st.E(2,0) * sv[0] + st.E(2,1) * sv[1] + st.E(2,2) * sv[2],
      st.E(0,0) * v_rxw[0] + st.E(0,1) * v_rxw[1] + st.E(0,2) * v_rxw[2],
      st.E(1,0) * v_rxw[0] + st.E(1,1) * v_rxw[1] + st.E(1,2) * v_rxw[2],
      st.E(2,0) * v_rxw[0] + st.E(2,1) * v_rxw[1] + st.E(2,2) * v_rxw[2]
      );
}

inline void applySTSV(
    unsigned ndirs,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    MatrixNd const & sv_dirs,
    SpatialVector & res,
    std::vector<SpatialVector> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.cols());
  assert(ndirs <= res_dirs.size());

  Vector3d v_rxw (
      sv[3] - st.r[1] * sv[2] + st.r[2] * sv[1],
      sv[4] - st.r[2] * sv[0] + st.r[0] * sv[2],
      sv[5] - st.r[0] * sv[1] + st.r[1] * sv[0]
      );

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d w_dir = sv_dirs.col(idir).segment<3>(0);
    res_dirs[idir].segment<3>(0) =
        st.E * w_dir
        + st_dirs[idir].E * sv.segment<3>(0);
    res_dirs[idir].segment<3>(3) =
        st.E * (sv_dirs.col(idir).segment<3>(3)
                - st.r.cross(w_dir)
                - st_dirs[idir].r.cross(sv.segment<3>(0)))
        + st_dirs[idir].E * v_rxw;
  }

  res = SpatialVector(
      st.E(0,0) * sv[0] + st.E(0,1) * sv[1] + st.E(0,2) * sv[2],
      st.E(1,0) * sv[0] + st.E(1,1) * sv[1] + st.E(1,2) * sv[2],
      st.E(2,0) * sv[0] + st.E(2,1) * sv[1] + st.E(2,2) * sv[2],
      st.E(0,0) * v_rxw[0] + st.E(0,1) * v_rxw[1] + st.E(0,2) * v_rxw[2],
      st.E(1,0) * v_rxw[0] + st.E(1,1) * v_rxw[1] + st.E(1,2) * v_rxw[2],
      st.E(2,0) * v_rxw[0] + st.E(2,1) * v_rxw[1] + st.E(2,2) * v_rxw[2]
      );
}

inline void applyTransposeSTSV (
    unsigned ndirs,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    std::vector<SpatialVector> const & sv_dirs,
    SpatialVector & res,
    std::vector<SpatialVector> & res_dirs
) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.size());
  assert(ndirs <= res_dirs.size());

  Vector3d E_T_f = st.E.transpose() * sv.segment<3>(3);
  Vector3d E_T_n = st.E.transpose() * sv.segment<3>(0);
  Vector3d top = E_T_n + st.r.cross(E_T_f);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d E_T_f_dir = // if &res == &sv
        st.E.transpose() * sv_dirs[idir].segment<3>(3)
        + st_dirs[idir].E.transpose() * sv.segment<3>(3);

    res_dirs[idir].segment<3>(0) =
        st.E.transpose() * sv_dirs[idir].segment<3>(0)
        + st_dirs[idir].E.transpose() * sv.segment<3>(0)
        + st_dirs[idir].r.cross(E_T_f)
        + st.r.cross(E_T_f_dir);

    res_dirs[idir].segment<3>(3) = E_T_f_dir;
  }

  res.segment<3>(0) = top;
  res.segment<3>(3) = E_T_f;
}

inline void addApplyTransposeSTSV (
    unsigned ndirs,
    double dscal,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    MatrixNd const & sv_dirs,
    SpatialVector & res,
    MatrixNd & res_dirs
) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.cols());
  assert(ndirs <= res_dirs.cols());

  Vector3d E_T_f = st.E.transpose() * sv.segment<3>(3);
  Vector3d E_T_n = st.E.transpose() * sv.segment<3>(0);
  Vector3d top = E_T_n + st.r.cross(E_T_f);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d E_T_f_dir = // if &res == &sv
        st.E.transpose() * sv_dirs.col(idir).segment<3>(3)
        + st_dirs[idir].E.transpose() * sv.segment<3>(3);

    res_dirs.col(idir).segment<3>(0) +=
        dscal * (
        st.E.transpose() * sv_dirs.col(idir).segment<3>(0)
        + st_dirs[idir].E.transpose() * sv.segment<3>(0)
        + st_dirs[idir].r.cross(E_T_f)
        + st.r.cross(E_T_f_dir));

    res_dirs.col(idir).segment<3>(3) += dscal * E_T_f_dir;
  }

  res.segment<3>(0) += dscal * top;
  res.segment<3>(3) += dscal * E_T_f;
}

inline void addApplyTransposeSTSV (
    unsigned ndirs,
    double dscal,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    std::vector<SpatialVector> const & sv_dirs,
    SpatialVector & res,
    std::vector<SpatialVector> & res_dirs
) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.size());
  assert(ndirs <= res_dirs.size());

  Vector3d E_T_f = st.E.transpose() * sv.segment<3>(3);
  Vector3d E_T_n = st.E.transpose() * sv.segment<3>(0);
  Vector3d top = E_T_n + st.r.cross(E_T_f);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d E_T_f_dir = // if &res == &sv
        st.E.transpose() * sv_dirs[idir].segment<3>(3)
        + st_dirs[idir].E.transpose() * sv.segment<3>(3);

    res_dirs[idir].segment<3>(0) +=
        dscal * (
        st.E.transpose() * sv_dirs[idir].segment<3>(0)
        + st_dirs[idir].E.transpose() * sv.segment<3>(0)
        + st_dirs[idir].r.cross(E_T_f)
        + st.r.cross(E_T_f_dir));

    res_dirs[idir].segment<3>(3) += dscal * E_T_f_dir;
  }

  res.segment<3>(0) += dscal * top;
  res.segment<3>(3) += dscal * E_T_f;
}

inline void addApplyAdjointSTSV (
    unsigned ndirs,
    double dscal,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    std::vector<SpatialVector> const & sv_dirs,
    SpatialVector & res,
    std::vector<SpatialVector> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.size());
  assert(ndirs <= res_dirs.size());

  Vector3d n_rxf = sv.segment<3>(0) - st.r.cross(sv.segment<3>(3));
  Vector3d En_rxf = st.E * n_rxf;
  Vector3d Ef = st.E * sv.segment<3>(3);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d En_rxf_dir = // if &res == &sv
        st_dirs[idir].E * n_rxf
        + st.E * (
          sv_dirs[idir].segment<3>(0)
          - st.r.cross(sv_dirs[idir].segment<3>(3))
          - st_dirs[idir].r.cross(sv.segment<3>(3))
          );

    res_dirs[idir].segment<3>(3) +=
        dscal * (
        st_dirs[idir].E * sv.segment<3>(3)
        + st.E * sv_dirs[idir].segment<3>(3));

    res_dirs[idir].segment<3>(0) += dscal * En_rxf_dir;
  }

  res.segment<3>(0) += dscal * En_rxf;
  res.segment<3>(3) += dscal * Ef;
}

inline void addApplyAdjointSTSV (
    unsigned ndirs,
    double dscal,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    MatrixNd const & sv_dirs,
    SpatialVector & res,
    std::vector<SpatialVector> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.cols());
  assert(ndirs <= res_dirs.size());

  Vector3d n_rxf = sv.segment<3>(0) - st.r.cross(sv.segment<3>(3));
  Vector3d En_rxf = st.E * n_rxf;
  Vector3d Ef = st.E * sv.segment<3>(3);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d En_rxf_dir = // if &res == &sv
        st_dirs[idir].E * n_rxf
        + st.E * (
          sv_dirs.col(idir).segment<3>(0)
          - st.r.cross(sv_dirs.col(idir).segment<3>(3))
          - st_dirs[idir].r.cross(sv.segment<3>(3))
          );

    res_dirs[idir].segment<3>(3) +=
        dscal * (
        st_dirs[idir].E * sv.segment<3>(3)
        + st.E * sv_dirs.col(idir).segment<3>(3));

    res_dirs[idir].segment<3>(0) += dscal * En_rxf_dir;
  }

  res.segment<3>(0) += dscal * En_rxf;
  res.segment<3>(3) += dscal * Ef;
}

inline void applyAdjointSTSV (
    unsigned ndirs,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    std::vector<SpatialVector> const & sv_dirs,
    SpatialVector & res,
    std::vector<SpatialVector> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.size());
  assert(ndirs <= res_dirs.size());

  Vector3d n_rxf = sv.segment<3>(0) - st.r.cross(sv.segment<3>(3));
  Vector3d En_rxf = st.E * n_rxf;
  Vector3d Ef = st.E * sv.segment<3>(3);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d En_rxf_dir = // if &res == &sv
        st_dirs[idir].E * n_rxf
        + st.E * (
          sv_dirs[idir].segment<3>(0)
          - st.r.cross(sv_dirs[idir].segment<3>(3))
          - st_dirs[idir].r.cross(sv.segment<3>(3))
          );

    res_dirs[idir].segment<3>(3) =
        st_dirs[idir].E * sv.segment<3>(3)
        + st.E * sv_dirs[idir].segment<3>(3);

    res_dirs[idir].segment<3>(0) = En_rxf_dir;
  }

  res.segment<3>(0) = En_rxf;
  res.segment<3>(3) = Ef;
}

inline void applyAdjointSTSV (
  unsigned ndirs,
  SpatialTransform const & st,
  std::vector<SpatialTransform> const & st_dirs,
  SpatialVector const & sv,
  MatrixNd const & sv_dirs,
  SpatialVector & res,
  MatrixNd & res_dirs
) {
  assert(sv_dirs.rows() == 6);
  assert(res_dirs.rows() == 6);
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.cols());
  assert(ndirs <= res_dirs.cols());

  Vector3d n_rxf = sv.segment<3>(0) - st.r.cross(sv.segment<3>(3));
  Vector3d En_rxf = st.E * n_rxf;
  Vector3d Ef = st.E * sv.segment<3>(3);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d En_rxf_dir = // if &res == &sv
        st_dirs[idir].E * n_rxf
        + st.E * (
          sv_dirs.col(idir).segment<3>(0)
          - st.r.cross(sv_dirs.col(idir).segment<3>(3))
          - st_dirs[idir].r.cross(sv.segment<3>(3))
          );

    res_dirs.col(idir).segment<3>(3) =
        st_dirs[idir].E * sv.segment<3>(3)
        + st.E * sv_dirs.col(idir).segment<3>(3);

    res_dirs.col(idir).segment<3>(0) = En_rxf_dir;
  }

  res.segment<3>(0) = En_rxf;
  res.segment<3>(3) = Ef;
}

inline void applyAdjointSTSV (
  unsigned ndirs,
  SpatialTransform const & st,
  std::vector<SpatialTransform> const & st_dirs,
  SpatialVector const & sv,
  SpatialVector & res,
  MatrixNd & res_dirs
) {
  assert(res_dirs.rows() == 6);
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= res_dirs.cols());

  Vector3d n_rxf = sv.segment<3>(0) - st.r.cross(sv.segment<3>(3));
  Vector3d En_rxf = st.E * n_rxf;
  Vector3d Ef = st.E * sv.segment<3>(3);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Vector3d En_rxf_dir = // if &res == &sv
        st_dirs[idir].E * n_rxf
        + st.E * (- st_dirs[idir].r.cross(sv.segment<3>(3)));

    res_dirs.col(idir).segment<3>(3) = st_dirs[idir].E * sv.segment<3>(3);
    res_dirs.col(idir).segment<3>(0) = En_rxf_dir;
  }

  res.segment<3>(0) = En_rxf;
  res.segment<3>(3) = Ef;
}

inline void applyTransposeSTSI (
    unsigned ndirs,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialRigidBodyInertia const & si,
    std::vector<SpatialRigidBodyInertia> const & si_dirs,
    SpatialRigidBodyInertia & res,
    std::vector<SpatialRigidBodyInertia> & res_dirs) {

  MatrixNd E_T_mr_dirs(3, ndirs);
  Vector3d E_T_mr;

  // derivative code
  for (unsigned i = 0; i < ndirs; i++) {
    E_T_mr_dirs.col(i) =
        st_dirs[i].E.transpose() * si.h + st.E.transpose() * si_dirs[i].h
        + si.m * st_dirs[i].r + si_dirs[i].m * st.r;
  }
  // nominal code
  E_T_mr = st.E.transpose() * si.h + si.m * st.r;

  Matrix3d Ixyz = si.getInertia();
  Matrix3d res_Ixyz_1 = st.E.transpose() * Ixyz * st.E;

  // derivative code
  for (unsigned i = 0; i < ndirs; i++) {
    res_dirs[i] = SpatialRigidBodyInertia (
        si_dirs[i].m,
        E_T_mr_dirs.col(i),
        st.E.transpose() * Ixyz * st_dirs[i].E
        + st.E.transpose() * si_dirs[i].getInertia() * st.E
        + st_dirs[i].E.transpose() * Ixyz * st.E

        - VectorCrossMatrix(st.r) * VectorCrossMatrix (st.E.transpose() * si_dirs[i].h)
        - VectorCrossMatrix(st.r) * VectorCrossMatrix (st_dirs[i].E.transpose() * si.h)
        - VectorCrossMatrix(st_dirs[i].r) * VectorCrossMatrix (st.E.transpose() * si.h)

        - VectorCrossMatrix (E_T_mr) * VectorCrossMatrix (st_dirs[i].r)
        - VectorCrossMatrix (E_T_mr_dirs.col(i)) * VectorCrossMatrix (st.r)
        );
  }

  // nominal code
  res = Math::SpatialRigidBodyInertia (
      si.m,
      E_T_mr,
      res_Ixyz_1
      - VectorCrossMatrix (st.r) * VectorCrossMatrix (st.E.transpose() * si.h)
      - VectorCrossMatrix (E_T_mr) * VectorCrossMatrix (st.r)
      );
}

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------


#endif
