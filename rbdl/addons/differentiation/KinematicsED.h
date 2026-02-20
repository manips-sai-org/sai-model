#ifndef RBDL_KINEMATICS_ED_h
#define RBDL_KINEMATICS_ED_h

#include "ModelED.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace ED {
// -----------------------------------------------------------------------------

/*
RBDL_DLLAPI Math::Vector3d CalcBodyToBaseCoordinates (
    Model &model,
    EDModel &ed_model,
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    unsigned int body_id,
    const Math::Vector3d &point_body_coordinates,
    Math::MatrixNd & ad_body_to_base_coordinates,
    bool update_kinematics = true);

RBDL_DLLAPI Math::Vector3d CalcBaseToBodyCoordinates (
    Model & model,
    EDModel & ed_model,
    Math::VectorNd const & q,
    Math::MatrixNd const & q_dirs,
    unsigned int body_id,
    Math::Vector3d const & base_point_position,
    Math::MatrixNd const & base_point_position_dirs,
    Math::MatrixNd & ad_base_to_body_coordinates,
    bool update_kinematics = true);
*/

RBDL_DLLAPI void UpdateKinematicsCustom (
    Model &model,
    EDModel & ed_model,
    Math::VectorNd const * q,
    Math::MatrixNd const * q_dirs,
    Math::VectorNd const * qdot,
    Math::MatrixNd const * qdot_dirs,
    Math::VectorNd const * qddot,
    Math::MatrixNd const * qddot_dirs);

RBDL_DLLAPI void UpdateKinematics (
    Model & model,
    EDModel & ed_model,
    Math::VectorNd const & q,
    Math::MatrixNd const & q_dirs,
    Math::VectorNd const & qdot,
    Math::MatrixNd const & qdot_dirs,
    Math::VectorNd const & qddot,
    Math::MatrixNd const & qddot_dirs);

RBDL_DLLAPI Math::Matrix3d CalcBodyWorldOrientation (
    Model & model,
    EDModel & ed_model,
    Math::VectorNd const & q,
    Math::MatrixNd const & q_dirs,
    const unsigned int body_id,
    std::vector<Math::Matrix3d> & ad_derivative,
    bool update_kinematics = true);

/*
RBDL_DLLAPI Math::Vector3d CalcPointVelocity (
    Model & model,
    EDModel & ed_model,
    Math::VectorNd const & q,
    Math::MatrixNd const & q_dirs,
    Math::VectorNd const & qdot,
    Math::MatrixNd const & qdot_dirs,
    unsigned int body_id,
    Math::Vector3d const & point_position,
    Math::MatrixNd & ad_point_velocity,
    bool update_kinematics = true);

RBDL_DLLAPI Math::SpatialVector CalcPointVelocity6D (
    Model &model,
    EDModel &ed_model,
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    const Math::VectorNd &qdot,
    const Math::MatrixNd &qdot_dirs,
    unsigned int body_id,
    const Math::Vector3d &point_position,
    std::vector<Math::SpatialVector> &pv6d_dirs,
    bool update_kinematics = true);
*/

RBDL_DLLAPI Math::Vector3d CalcPointAcceleration (
    Model & model,
    EDModel & ed_model,
    Math::VectorNd const & q,
    Math::MatrixNd const & q_dirs,
    Math::VectorNd const & qdot,
    Math::MatrixNd const & qdot_dirs,
    Math::VectorNd const & qddot,
    Math::MatrixNd const & qddot_dirs,
    unsigned int body_id,
    Math::Vector3d const & point_position,
    Math::MatrixNd & ad_derivative,
    bool update_kinematics = true);

RBDL_DLLAPI void CalcPointJacobian (
    Model &model,
    EDModel &ed_model,
    Math::VectorNd const &Q,
    Math::MatrixNd const &Q_dirs,
    unsigned int body_id,
    Math::Vector3d const &point_position,
    Math::MatrixNd &G,
    std::vector<Math::MatrixNd> &G_dirs,
    bool update_kinematics = true
    );

/*
RBDL_DLLAPI void CalcPointJacobian6D (
    Model &model,
    EDModel &ed_model,
    Math::VectorNd const &q,
    Math::MatrixNd const &q_dirs,
    unsigned body_id,
    Math::Vector3d const &point_position,
    Math::MatrixNd &G,
    std::vector<Math::MatrixNd> &ad_G,
    bool update_kinematics);
*/

// -----------------------------------------------------------------------------
} // namespace ED
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif
