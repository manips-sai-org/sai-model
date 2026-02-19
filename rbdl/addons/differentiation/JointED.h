#ifndef RBDL_JOINT_ED_h
#define RBDL_JOINT_ED_h

#include <rbdl/Model.h>

#include "ModelED.h"
#include "ModelAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace ED {
// -----------------------------------------------------------------------------

RBDL_DLLAPI void jcalc (
  RigidBodyDynamics::Model &model,
  EDModel &ed_model,
  unsigned int joint_id,
  const Math::VectorNd &q,
  const Math::MatrixNd &q_dirs,
  const Math::VectorNd &qdot,
  const Math::MatrixNd &qdot_dirs
);

RBDL_DLLAPI void jcalc_XJ (
  RigidBodyDynamics::Model &model,
  EDModel &ed_model,
  const unsigned int &joint_id,
  const unsigned int &ndir,
  const Math::VectorNd &q,
  const Math::MatrixNd &q_dirs,
  Math::SpatialTransform &X,
  std::vector<Math::SpatialTransform> &X_dir
);

// -----------------------------------------------------------------------------
} // namespace ED
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif
