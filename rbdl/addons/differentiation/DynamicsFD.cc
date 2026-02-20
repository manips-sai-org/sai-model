/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */
#include "DynamicsFD.h"
#include "FdModelEntry.h"

using namespace RigidBodyDynamics::Math;

using std::vector;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

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
         && "q_dirs, qdot_dirs, tau_dirs, fd_qddot have different dimensions");

  double h = sqrt(1e-16);
  unsigned int ndirs = q_dirs.cols();

  ForwardDynamics(model, q, qdot, tau, qddot,f_ext);

  VectorNd hd_qddot(qddot);
  VectorNd q_dir(q);
  VectorNd qdot_dir(qdot);
  VectorNd tau_dir(tau);

  for(unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }
    VectorNd qh  = q + h * q_dirs.col(idir);
    VectorNd qdh = qdot + h * qdot_dirs.col(idir);
    VectorNd tauh = tau + h * tau_dirs.col(idir);
    vector<SpatialVector> f_exth_;
    if (f_ext) {
      f_exth_.resize(f_ext->size());
      for (unsigned i = 0; i < f_ext->size(); i++) {
        f_exth_[i] = (*f_ext)[i] + h * (*f_ext_dirs)[i][idir];
      }
    }
    vector<SpatialVector> const * f_exth = f_ext ? &f_exth_ : NULL;
    ForwardDynamics(*modelh, qh, qdh, tauh, hd_qddot, f_exth);
    fd_qddot.col(idir) = (hd_qddot - qddot) / h;
    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
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

  double h = sqrt(1e-16);

  unsigned int ndirs   = q_dirs.cols();

  InverseDynamics(model, q, qdot, qddot, tau, f_ext);

  VectorNd tau_h (tau);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }

    vector<SpatialVector> f_exth_;
    if (f_ext) {
      f_exth_.resize(f_ext->size());
      for (unsigned i = 0; i < f_ext->size(); i++) {
        f_exth_[i] = (*f_ext)[i] + h * (*f_ext_dirs)[i][idir];
      }
    }

    InverseDynamics(
      *modelh, q + h*q_dirs.col(idir),
      qdot + h*qdot_dirs.col(idir),
      qddot + h*qddot_dirs.col(idir),
      tau_h,
      f_ext ? &f_exth_ : NULL
    );

    fd_tau.col(idir) = (tau_h - tau) / h;

    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
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
  unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == fd_tau.cols());

  NonlinearEffects (model, q, qdot, tau);
  double h = 1.0e-8;
  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }
    VectorNd q_dir  = q_dirs.col(idir);
    VectorNd qd_dir = qdot_dirs.col(idir);
    VectorNd hd_tau(tau.rows());
    NonlinearEffects (*modelh, q + h * q_dir, qdot + h * qd_dir, hd_tau);
    fd_tau.col(idir) = (hd_tau - tau) / h;
    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
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

  // value of perturbation
  double h = 1.0e-8;

  // compute value at current configuration q as nominal value
  CompositeRigidBodyAlgorithm(model, q, H, true);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    // hack to compute finite differences of model quantities
    Model *modelh = &model;
    if (fd_model) {
      modelh = new Model(model);
    }
    VectorNd q_dir = q_dirs.col(idir);

    // compute perturbations along direction
    CompositeRigidBodyAlgorithm(*modelh, q+h*q_dir, fd_H[idir], true);

    // compute directional derivatives using forward first order differences
    fd_H[idir] = (fd_H[idir] - H) / h;

    // hack to compute finite differences of model quantities
    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
}

// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

