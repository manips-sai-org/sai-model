/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_MODEL_AD_H
#define RBDL_MODEL_AD_H

#include <assert.h>
#include <iostream>
#include <limits>
#include <cstring>

#include "rbdl/rbdl.h"
#include "rbdl/rbdl_math.h"

//using namespace std;
//using namespace RigidBodyDynamics;
//using namespace RigidBodyDynamics::Math;

// // -----------------------------------------------------------------------------
// namespace RigidBodyDynamics {
// // -----------------------------------------------------------------------------
// namespace AD {
// -----------------------------------------------------------------------------

struct ADModel {

    unsigned int ndirs;

    // derivative values
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialTransform> > X_lambda;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialTransform> > X_base;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialTransform> > X_J;

    std::vector<RigidBodyDynamics::Math::SpatialTransform> X_0;

    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > v_J;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > v;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > c_J;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > a;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > c;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > f;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > h;

    RigidBodyDynamics::Math::SpatialRigidBodyInertia I0;
    std::vector<RigidBodyDynamics::Math::SpatialRigidBodyInertia> I_temp;

    std::vector<std::vector<RigidBodyDynamics::Math::SpatialRigidBodyInertia> > Ic;

    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > hc;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > F;

    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > pA;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialVector> > U;
    std::vector<std::vector<RigidBodyDynamics::Math::SpatialMatrix> > IA;

    //std::vector<std::vector<double> > u;
    //std::vector<std::vector<double> > d;
    RigidBodyDynamics::Math::MatrixNd u;
    RigidBodyDynamics::Math::MatrixNd d;

    ADModel ();
    ADModel (RigidBodyDynamics::Model& model);

    void resize_directions (const unsigned & requested_ndirs);

};

// -----------------------------------------------------------------------------
// } /* AD */
// -----------------------------------------------------------------------------
// } /* RigidBodyDynamics */
// -----------------------------------------------------------------------------

/* RBDL_MODEL_AD_H */
#endif
