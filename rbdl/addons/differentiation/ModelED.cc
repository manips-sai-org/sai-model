/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include "ModelED.h"

using namespace RigidBodyDynamics::Math;

// std::vectors containing any objectst that have Eigen matrices or vectors
// as members need to have a special allocater. This can be achieved with
// the following macro.

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
// namespace ED {
// -----------------------------------------------------------------------------

EDModel::EDModel (RigidBodyDynamics::Model const & model) {
    ndirs = 10 * model.dof_count;

    bodyidx2s1idx = std::vector<uint8_t>(model.mBodies.size(), 0);
    for (unsigned i = 0; i < model.mBodies.size(); i++) {
      for (unsigned j = 0; j < 6; j++) {
        if (model.S[i](j) == 1.) {
          bodyidx2s1idx[i] = j;
          break;
        }
      }
    }

    // NOTE: old initialization values
    // std::vector<SpatialRigidBodyInertia> SRBI(ndirs);
    std::vector<SpatialRigidBodyInertia> SRBI(ndirs);
    std::vector<SpatialMatrix> X(ndirs, SpatialMatrix::Zero());
    std::vector<SpatialTransform> T(ndirs, SpatialTransform::Zero());
    // std::vector<SpatialTransform> Tdot(ndirs, SpatialTransform::Zero());
    std::vector<SpatialVector> vec(ndirs, SpatialVector::Zero());
    SpatialDirection vec_dir = MatrixNd::Zero(6, ndirs);

    ad_E.resize(ndirs, Matrix3dZero);
    ad_p_X_i.resize(ndirs, SpatialMatrix::Zero());
    ad_p_v_i.resize(6, ndirs);
    ad_a_dash.resize(3, ndirs);
    ad_p_a_i.resize(6, ndirs);

    //std::vector<double> double(ndirs, 0.0);

    // inertias
    Ic.resize(model.mBodies.size(), SRBI);

    // std::cout << "init m bodies: " << model.mBodies.size() << "\n";
    // std::cout << "init Ic size: " << Ic.size() << "\n";

    // spatial matrices
    IA.resize(model.mBodies.size(), X);

    // spatial transforms
    X_lambda.resize(model.mBodies.size(), T);
    X_base.resize(model.mBodies.size(), T);
    X_J.resize(model.mBodies.size(), T);
    X_0.resize(ndirs, SpatialTransform::Zero());

    // spatial vectors
    v_J.resize(model.mBodies.size(), vec_dir);
    v.resize(model.mBodies.size(), vec_dir);
    v_q.resize(model.mBodies.size(), vec_dir);
    v_qdot.resize(model.mBodies.size(), vec_dir);
    v_qddot.resize(model.mBodies.size(), vec_dir);
    c_J.resize(model.mBodies.size(), vec_dir);
    a.resize(model.mBodies.size(), vec_dir);
    a_q.resize(model.mBodies.size(), vec_dir);
    a_qdot.resize(model.mBodies.size(), vec_dir);
    a_qddot.resize(model.mBodies.size(), vec_dir);

    SpatialVector temp = SpatialVector::Zero();
    h.resize(model.mBodies.size(), temp);
    h_q.resize(model.mBodies.size(), vec_dir);
    h_qdot.resize(model.mBodies.size(), vec_dir);
    h_qddot.resize(model.mBodies.size(), vec_dir);

    hc.resize(model.mBodies.size(), vec_dir);
    F.resize(model.mBodies.size(), vec_dir);
    c.resize(model.mBodies.size(), vec_dir);
    c_q.resize(model.mBodies.size(), vec_dir);
    c_qdot.resize(model.mBodies.size(), vec_dir);
    c_qddot.resize(model.mBodies.size(), vec_dir);
    f.resize(model.mBodies.size(), vec_dir);
    f_q.resize(model.mBodies.size(), vec_dir);
    f_qdot.resize(model.mBodies.size(), vec_dir);
    f_qddot.resize(model.mBodies.size(), vec_dir);

    pA.resize(model.mBodies.size(), vec_dir);
    U.resize(model.mBodies.size(), vec_dir);

    // spatial vectors
    u.resize(model.u.rows(), ndirs);
    u.setZero();
    d.resize(model.d.rows(), ndirs);
    d.setZero();

    Iv = vec_dir;
}

void EDModel::resize_directions (const unsigned int &requested_ndirs)
{
  if (ndirs < requested_ndirs) {
      ndirs = requested_ndirs;

    // inertias
    for (unsigned int i = 0; i < Ic.size(); ++i) {
        Ic[i].resize(ndirs);
    }

    // spatial matrices
    for (unsigned int i = 0; i < IA.size(); ++i) {
        IA[i].resize(ndirs, SpatialMatrix::Zero());
    }

    // spatial transforms
    for (unsigned int i = 0; i < X_lambda.size(); ++i) {
        X_lambda[i].resize(ndirs, SpatialTransform::Zero());
    }

    for (unsigned int i = 0; i < X_base.size(); ++i) {
        X_base[i].resize(ndirs, SpatialTransform::Zero());
    }

    for (unsigned int i = 0; i < X_J.size(); ++i) {
        X_J[i].resize(ndirs, SpatialTransform::Zero());
    }

    X_0.resize(ndirs, SpatialTransform::Zero());

    // spatial vectors
    for (unsigned int i = 0; i < v_J.size(); ++i) {
        v_J[i].resize(6, ndirs);
    }

    for (unsigned int i = 0; i < v_J.size(); ++i) {
        v[i].resize(6, ndirs);
    }

    for (unsigned int i = 0; i < c_J.size(); ++i) {
        c_J[i].resize(6, ndirs);
    }

    for (unsigned int i = 0; i < c_J.size(); ++i) {
        a[i].resize(6, ndirs);
    }

    for (unsigned int i = 0; i < c_J.size(); ++i) {
        c[i].resize(6, ndirs);
    }

    for (unsigned int i = 0; i < c_J.size(); ++i) {
        f[i].resize(6, ndirs);
    }

    for (unsigned int i = 0; i < F.size(); ++i) {
        F[i].resize(6, ndirs);
    }

    for(unsigned int i = 0; i < pA.size(); ++i) {
        pA[i].resize(6, ndirs);
    }

    for(unsigned int i = 0; i < U.size(); ++i) {
        U[i].resize(6, ndirs);
    }

    for (unsigned int i = 0; i < hc.size(); ++i) {
        hc[i].resize(6, ndirs);
    }

    // other quantities
    u.resize(u.rows(), ndirs);
    u.setZero();
    d.resize(d.rows(), ndirs);
    d.setZero();
    Iv.resize(6, ndirs);
    ad_E.resize(ndirs, Matrix3dZero);
    ad_p_X_i.resize(ndirs, SpatialMatrix::Zero());
    ad_p_v_i.resize(6, ndirs);
    ad_a_dash.resize(3, ndirs);
    ad_p_a_i.resize(6, ndirs);
  }
}

EDModel::~EDModel () {};

// -----------------------------------------------------------------------------
// } // ED
// -----------------------------------------------------------------------------
} // RigidBodyDynamics
// -----------------------------------------------------------------------------
