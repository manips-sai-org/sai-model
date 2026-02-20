/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_DYNAMICS_FD_H
#define RBDL_DYNAMICS_FD_H

#include <iostream>
#include <limits>
#include <assert.h>
#include <string.h>

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Joint.h"
#include "rbdl/Body.h"
#include "rbdl/Dynamics.h"
#include "rbdl/Kinematics.h"
#include "ModelAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void ForwardDynamics (
    Model &model,
    ADModel *fd_model, // NULL means execution without fd_model update
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    const Math::VectorNd &qdot,
    const Math::MatrixNd &qdot_dirs,
    const Math::VectorNd &tau,
    const Math::MatrixNd &tau_dirs,
    Math::VectorNd &qddot,
    Math::MatrixNd &fd_qddot,
    std::vector<Math::SpatialVector> const *f_ext = NULL,
    std::vector<std::vector<Math::SpatialVector> > const *f_ext_dirs = NULL);

RBDL_DLLAPI
void InverseDynamics (
    Model &model,
    ADModel *fd_model, // NULL means execution without fd_model update
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    const Math::VectorNd &qdot,
    const Math::MatrixNd &qdot_dirs,
    const Math::VectorNd &qddot,
    const Math::MatrixNd &qddot_dirs,
    Math::VectorNd &tau,
    Math::MatrixNd &fd_tau,
    std::vector<Math::SpatialVector> const *f_ext = NULL,
    std::vector<std::vector<Math::SpatialVector> > const *f_ext_dirs = NULL);

RBDL_DLLAPI
void NonlinearEffects (
    Model &model,
    ADModel *fd_model, // NULL means execution without fd_model update
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    const Math::VectorNd &qdot,
    const Math::MatrixNd &qdot_dirs,
    Math::VectorNd &tau,
    Math::MatrixNd &fd_tau);

RBDL_DLLAPI
void CompositeRigidBodyAlgorithm (
    Model &model,
    ADModel *fd_model,
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    Math::MatrixNd &H,
    std::vector<Math::MatrixNd> &fd_H);

// -----------------------------------------------------------------------------
} /* FD */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------

/* RBDL_DYNAMICS_FD_H */
#endif
