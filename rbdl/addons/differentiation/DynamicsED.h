/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_DYNAMICS_ED_H
#define RBDL_DYNAMICS_ED_H


// #include "ModelAD.h"
#include "ModelED.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace ED {
// -----------------------------------------------------------------------------

RBDL_DLLAPI void InverseDynamics(
  Model &model,
  EDModel &ad_model,
  Math::VectorNd const &q,
  Math::MatrixNd const &q_dirs,
  Math::VectorNd const &qdot,
  Math::MatrixNd const &qdot_dirs,
  Math::VectorNd const &qddot,
  Math::MatrixNd const &qddot_dirs,
  Math::VectorNd &tau,
  Math::MatrixNd &ad_tau,
  std::vector<Math::SpatialVector> const *f_ext = nullptr,
  std::vector<std::vector<Math::SpatialVector> > const *f_ext_dirs = nullptr
);


RBDL_DLLAPI
void NonlinearEffects (
    Model & model,
    EDModel & ed_model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot,
    const Math::MatrixNd & qdot_dirs,
    Math::VectorNd & tau,
    Math::MatrixNd & ed_tau
    );


RBDL_DLLAPI
void CompositeRigidBodyAlgorithm (
  Model &model,
  EDModel &ed_model,
  Math::VectorNd const & q,
  Math::MatrixNd const & q_dirs,
  Math::MatrixNd & H,
  std::vector<Math::MatrixNd> & H_dirs,
  bool update_kinematics = true
);



// -----------------------------------------------------------------------------
} /* ED */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------

/* RBDL_DYNAMICS_ED_H */
#endif
