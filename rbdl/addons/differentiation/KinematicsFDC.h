#ifndef RBDL_KINEMATICS_FDC_h
#define RBDL_KINEMATICS_FDC_h

#include <rbdl/Model.h>
#include <rbdl/rbdl_math.h>

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FDC {
// -----------------------------------------------------------------------------

RBDL_DLLAPI Math::Vector3d CalcBodyToBaseCoordinatesSingleFunc (
		Model &model,
		const Math::VectorNd &q,
		const Math::MatrixNd &q_dirs,
		unsigned int body_id,
		const Math::Vector3d &point_body_coordinates,
		Math::MatrixNd &out);

RBDL_DLLAPI Math::Vector3d CalcBodyToBaseCoordinates (
    Model & model,
    ADModel * fd_model,
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    unsigned int body_id,
    const Math::Vector3d &point_body_coordinates,
    Math::MatrixNd &fd_body_to_base_coordinates);

RBDL_DLLAPI Math::Vector3d CalcBaseToBodyCoordinates (
    Model & model,
    ADModel * fd_model,
    Math::VectorNd const & q,
    Math::MatrixNd const & q_dirs,
    unsigned body_id,
    Math::Vector3d const & base_point_position,
    Math::MatrixNd const & base_point_position_dirs,
    Math::MatrixNd & fd_base_to_body_coordinates);

RBDL_DLLAPI Math::Matrix3d CalcBodyWorldOrientation (
    Model & model,
    Math::VectorNd const & q,
    Math::MatrixNd const & q_dirs,
    unsigned int const body_id,
    std::vector<Math::Matrix3d> & fd_derivative);

RBDL_DLLAPI Math::Vector3d CalcPointVelocity (
    Model &model,
    Math::VectorNd const &q,
    Math::MatrixNd const &q_dirs,
    Math::VectorNd const &qdot,
    Math::MatrixNd const &qdot_dirs,
    unsigned int body_id,
    Math::Vector3d const &point_position,
    Math::MatrixNd &fd_point_velocity);

RBDL_DLLAPI Math::SpatialVector CalcPointVelocity6D (
    Model &model,
    ADModel *fd_model, // NULL means execution without fd_model update
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    const Math::VectorNd &qdot,
    const Math::MatrixNd &qdot_dirs,
    unsigned int body_id,
    const Math::Vector3d &point_position,
    std::vector<Math::SpatialVector> &pv6d_dirs);

RBDL_DLLAPI Math::Vector3d CalcPointAcceleration (
    Model & model,
    Math::VectorNd const & q,
    Math::MatrixNd const & q_dirs,
    Math::VectorNd const & qdot,
    Math::MatrixNd const & qdot_dirs,
    Math::VectorNd const & qddot,
    Math::MatrixNd const & qddot_dirs,
    unsigned int body_id,
    Math::Vector3d const & point_position,
    Math::MatrixNd & fd_derivative);

RBDL_DLLAPI void CalcPointJacobian (
    Model &model,
    ADModel *fd_model,
    Math::VectorNd const &q,
    Math::MatrixNd const &q_dirs,
    unsigned int body_id,
    Math::Vector3d const &point_position,
    Math::MatrixNd &G,
    std::vector<Math::MatrixNd> &G_dirs);

RBDL_DLLAPI void CalcPointJacobian6D (
    Model &model,
    ADModel *fd_model, // NULL means execution without fd_model update
    Math::VectorNd const &q,
    Math::MatrixNd const &q_dirs,
    unsigned body_id,
    Math::Vector3d const &point_position,
    Math::MatrixNd &G,
    std::vector<Math::MatrixNd> &fd_G);

RBDL_DLLAPI void UpdateKinematicsCustom (
    Model & model,
    ADModel *fd_model,
    Math::VectorNd const & q,
    Math::MatrixNd const & q_dirs,
    Math::VectorNd const & qd,
    Math::MatrixNd const & qd_dirs,
    Math::VectorNd const & qdd,
    Math::MatrixNd const & qdd_dirs);

// -----------------------------------------------------------------------------
} // namespace FDC
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif
