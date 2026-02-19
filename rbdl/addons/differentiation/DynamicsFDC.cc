/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */
#include "DynamicsFDC.h"
#include "ModelEntryFDC.h"

#include <cfloat>

using std::vector;
using namespace RigidBodyDynamics::Math;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FDC {
// -----------------------------------------------------------------------------

// define constant difference perturbation according to theory,
// i.e. EPS = eps^(1/3) = cubic_root(eps).
// Here, we us cmath cubic root implementation, i.e. cubic_root = cbrt.
// NOTE assumes that directions are normalized to one
const double EPS = cbrt(DBL_EPSILON);
const double EPSx2 = 2*EPS; // for convenience, we directly compute the denominator

RBDL_DLLAPI
void ForwardDynamics(
    Model &model,
    ADModel *fd_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs,
    const VectorNd &tau,
    const MatrixNd &tau_dirs,
    VectorNd &qddot,
    MatrixNd &fd_qddot,
    vector<SpatialVector> const *f_ext,
    vector<vector<SpatialVector>> const *f_ext_dirs
) {
  assert((f_ext == NULL) == (f_ext_dirs == NULL));

  assert(q_dirs.cols() == qdot_dirs.cols()
    && q_dirs.cols() == fd_qddot.cols()
    && tau_dirs.cols() == q_dirs.cols()
    && "q_dirs, qdot_dirs, tau_dirs, fd_qddot have different dimensions"
  );

  const unsigned int ndirs = q_dirs.cols();

  // temporary results
  VectorNd qddot_ph(qddot);
  VectorNd qddot_mh(qddot);

  // handle external forces
  vector<SpatialVector> f_exth;
  if (f_ext) {
    f_exth.resize(f_ext->size());
  }

  // compute central difference
  for(unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }

    // evaluate forward
    if (f_ext) {
      for (unsigned i = 0; i < f_ext->size(); i++) {
        f_exth[i] = (*f_ext)[i] + EPS * (*f_ext_dirs)[i][idir];
      }
    }

    ForwardDynamics(
      *modelh,
      q + EPS * q_dirs.col(idir),
      qdot + EPS * qdot_dirs.col(idir),
      tau + EPS * tau_dirs.col(idir),
      qddot_ph,
      f_ext ? &f_exth : NULL
    );

    // evaluate backward
    if (f_ext) {
      for (unsigned i = 0; i < f_ext->size(); i++) {
        f_exth[i] = (*f_ext)[i] - EPS * (*f_ext_dirs)[i][idir];
      }
    }

    ForwardDynamics(
      model,
      q - EPS * q_dirs.col(idir),
      qdot - EPS * qdot_dirs.col(idir),
      tau - EPS * tau_dirs.col(idir),
      qddot_mh,
      f_ext ? &f_exth : NULL
    );

    fd_qddot.col(idir) = (qddot_ph - qddot_mh) / EPSx2;
    if (fd_model) {
      computeFDCEntry(*modelh, model, EPS, idir, *fd_model);
      delete modelh;
    }
  }

  // nominal evaluation
  ForwardDynamics(model, q, qdot, tau, qddot,f_ext);
}



RBDL_DLLAPI
void InverseDynamics(
    Model &model,
    ADModel *fd_model,
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    const Math::VectorNd &qdot,
    const Math::MatrixNd &qdot_dirs,
    const Math::VectorNd &qddot,
    const Math::MatrixNd &qddot_dirs,
    Math::VectorNd &tau,
    Math::MatrixNd &fd_tau,
    vector<SpatialVector> const *f_ext,
    vector<vector<SpatialVector>> const *f_ext_dirs
) {
  assert(q_dirs.cols() == qdot_dirs.cols() &&
         q_dirs.cols() == qddot_dirs.cols() &&
         "q_dirs, qdot_dirs, qddot_dirs have different dimensions");

  assert( fd_tau.cols() == q_dirs.cols() &&
          fd_tau.rows() == tau.rows() &&
          "fd_tau and tau have different dimensions");


  const unsigned int ndirs   = q_dirs.cols();

  // handling external forces
  vector<SpatialVector> f_exth;
  if (f_ext) {
    f_exth.resize(f_ext->size());
  }

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }

    // forward perturbation
    if (f_ext) {
      for (unsigned i = 0; i < f_ext->size(); i++) {
        f_exth[i] = (*f_ext)[i] + EPS * (*f_ext_dirs)[i][idir];
      }
    }

    InverseDynamics(
      *modelh,
      q + EPS*q_dirs.col(idir),
      qdot + EPS*qdot_dirs.col(idir),
      qddot + EPS*qddot_dirs.col(idir),
      tau,
      f_ext ? &f_exth : NULL
    );
    fd_tau.col(idir) = tau;

    // backward evaluation
    if (f_ext) {
      for (unsigned i = 0; i < f_ext->size(); i++) {
        f_exth[i] = (*f_ext)[i] - EPS * (*f_ext_dirs)[i][idir];
      }
    }

    InverseDynamics(
      model,
      q - EPS*q_dirs.col(idir),
      qdot - EPS*qdot_dirs.col(idir),
      qddot - EPS*qddot_dirs.col(idir),
      tau,
      f_ext ? &f_exth : NULL
    );

    fd_tau.col(idir) = (fd_tau.col(idir) - tau) / EPSx2;

    if (fd_model) {
      computeFDCEntry(*modelh, model , EPS, idir, *fd_model);
      delete modelh;
    }
  }

  // nominal evaluation
  InverseDynamics(model, q, qdot, qddot, tau, f_ext);

  return;
}

RBDL_DLLAPI
void NonlinearEffects (
    Model &model,
    ADModel *fd_model,
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    const Math::VectorNd &qdot,
    const Math::MatrixNd &qdot_dirs,
    Math::VectorNd &tau,
    Math::MatrixNd &fd_tau
) {
  const unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == fd_tau.cols());

  VectorNd tau_ph(tau.rows());
  VectorNd tau_mh(tau.rows());

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }

    // forward perturbation
    NonlinearEffects (
      *modelh,
      q + EPS * q_dirs.col(idir),
      qdot + EPS * qdot_dirs.col(idir),
      tau_ph
    );

    // backward perturbation
    NonlinearEffects (
      model,
      q - EPS * q_dirs.col(idir),
      qdot - EPS * qdot_dirs.col(idir),
      tau_mh
    );

    fd_tau.col(idir) = (tau_ph - tau_mh) / EPSx2;

    if (fd_model) {
      computeFDCEntry(*modelh, model, EPS, idir, *fd_model);
      delete modelh;
    }
  }

  NonlinearEffects (model, q, qdot, tau);
}

RBDL_DLLAPI
void CompositeRigidBodyAlgorithm (
    Model &model,
    ADModel *fd_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    MatrixNd &H,
    vector<MatrixNd> &fd_H
) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == fd_H.size());

  for (unsigned idir = 0; idir < ndirs; idir++) {
    // hack to compute finite differences of model quantities
    Model *modelh = &model;
    if (fd_model) {
      modelh = new Model(model);
    }

    // forward perturbation
    CompositeRigidBodyAlgorithm(
      *modelh,
      q + EPS * q_dirs.col(idir),
      fd_H[idir],
      true
     );

    CompositeRigidBodyAlgorithm(
      model,
      q - EPS * q_dirs.col(idir),
      H,
      true
     );

    // compute directional derivatives using forward first order differences
    fd_H[idir] = (fd_H[idir] - H) / EPSx2;

    // hack to compute finite differences of model quantities
    if (fd_model) {
      computeFDCEntry(*modelh, model, EPS, idir, *fd_model);
      delete modelh;
    }

  }
  // compute value at current configuration q as nominal value
  CompositeRigidBodyAlgorithm(model, q, H, true);

}


// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

