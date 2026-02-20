#ifndef RBDL_UTILS_AD_h
#define RBDL_UTILS_AD_h

#include "ModelAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {  
// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI void CalcCenterOfMass (
        Model & model,
        ADModel & ad_model,
        Math::VectorNd const & q,
        Math::MatrixNd const & q_dirs,
        Math::VectorNd const & qdot,
        Math::MatrixNd const & qdot_dirs,
        double & mass,
        Math::Vector3d & com,
        Math::MatrixNd & ad_com,
        Math::Vector3d * com_velocity,
        Math::MatrixNd * ad_com_velocity,
        Math::Vector3d * angular_momentum,
        Math::MatrixNd * ad_angular_momentum,
        bool update_kinematics = true);

RBDL_DLLAPI double CalcPotentialEnergy (
        Model & model,
        ADModel & ad_model,
        Math::VectorNd const & q,
        Math::MatrixNd const & q_dirs,
        Math::MatrixNd & ad_pote,
        bool update_kinematics = true);

RBDL_DLLAPI double CalcKineticEnergy (
        Model & model,
        ADModel & ad_model,
        Math::VectorNd const & q,
        Math::MatrixNd const & q_dirs,
        Math::VectorNd const & qdot,
        Math::MatrixNd const & qdot_dirs,
        Math::MatrixNd & ad_kine,
        bool update_kinematics = true);

// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------



#endif
