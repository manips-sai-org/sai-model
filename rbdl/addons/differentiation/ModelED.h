/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_MODEL_ED_H
#define RBDL_MODEL_ED_H

#include <assert.h>
#include <iostream>
#include <limits>
#include <cstring>

#include "rbdl/rbdl.h"
#include "rbdl/rbdl_math.h"

#include "SpatialAlgebraOperatorsED.h"

//using namespace std;
//using namespace RigidBodyDynamics;
//using namespace RigidBodyDynamics::Math;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
// namespace ED {
// -----------------------------------------------------------------------------

typedef Eigen::Matrix<double, 6, Eigen::Dynamic> SpatialDirection;

struct EDModel {

    unsigned int ndirs;

    std::vector<uint8_t> bodyidx2s1idx;
    std::vector<RigidBodyDynamics::Math::Matrix3d> ad_E;
    std::vector<RigidBodyDynamics::Math::SpatialMatrix> ad_p_X_i;
    RigidBodyDynamics::Math::MatrixNd ad_p_v_i;
    RigidBodyDynamics::Math::MatrixNd ad_a_dash;
    RigidBodyDynamics::Math::MatrixNd ad_p_a_i;

    // inertias
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialRigidBodyInertia> > Ic;

    // spatial matrices
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialMatrix> > IA;

    // spatial transforms
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialTransform> > X_lambda;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialTransform> > X_base;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialTransform> > X_J;
    std::vector<RigidBodyDynamics::Math::SpatialTransform> X_0;


    // spatial vectors
    std::vector<SpatialDirection> v_J;
    std::vector<SpatialDirection> v;
    std::vector<SpatialDirection> v_q;
    std::vector<SpatialDirection> v_qdot;
    std::vector<SpatialDirection> v_qddot;
    std::vector<SpatialDirection> c_J;
    std::vector<SpatialDirection> a;
    std::vector<SpatialDirection> a_q;
    std::vector<SpatialDirection> a_qdot;
    std::vector<SpatialDirection> a_qddot;
    std::vector<SpatialDirection> c;
    std::vector<SpatialDirection> c_q;
    std::vector<SpatialDirection> c_qdot;
    std::vector<SpatialDirection> c_qddot;
    std::vector<SpatialDirection> f;
    std::vector<SpatialDirection> f_q;
    std::vector<SpatialDirection> f_qdot;
    std::vector<SpatialDirection> f_qddot;

    std::vector<RigidBodyDynamics::Math::SpatialVector> h;
    std::vector<SpatialDirection> h_q;
    std::vector<SpatialDirection> h_qdot;
    std::vector<SpatialDirection> h_qddot;

    std::vector<SpatialDirection> hc;
    std::vector<SpatialDirection> F;

    std::vector<SpatialDirection> pA;
    std::vector<SpatialDirection> U;

    // other quantities
    RigidBodyDynamics::Math::MatrixNd u;
    RigidBodyDynamics::Math::MatrixNd d;
    SpatialDirection Iv;

    EDModel () = delete;
    EDModel (RigidBodyDynamics::Model const& model);
    ~EDModel ();

    void resize_directions (const unsigned int &requested_ndirs);

};

// -----------------------------------------------------------------------------
// } // ED
// -----------------------------------------------------------------------------
} // RigidBodyDynamics
// -----------------------------------------------------------------------------

/* RBDL_MODEL_ED_H */
#endif
