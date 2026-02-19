#ifndef RBDL_CONTACTS_FD_H
#define RBDL_CONTACTS_FD_H

#include <rbdl/rbdl_math.h>
#include <rbdl/rbdl_mathutils.h>

#include "rbdl/Constraints.h"

#include "KinematicsAD.h"
#include "ConstraintsAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI void CalcConstrainedSystemVariables (
    Model &model,
    ADModel *fd_model, // NULL means execution without fd_model update
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    const Math::VectorNd &qdot,
    const Math::MatrixNd &qdot_dirs,
    const Math::VectorNd &tau,
    const Math::MatrixNd &tau_dirs,
    ConstraintSet   &cs,
    ADConstraintSet *fd_CS);

RBDL_DLLAPI
void CalcConstraintsJacobian (
    Model &model,
    ADModel *fd_model, // NULL means execution without fd_model update
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    ConstraintSet &CS,
    ADConstraintSet &fd_CS,
    Math::MatrixNd &G,
    std::vector<Math::MatrixNd> &G_dirs);

RBDL_DLLAPI
void ComputeConstraintImpulsesDirect (
    Model &model,
    ADModel *fd_model,// NULL means execution without fd_model update
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot_minus,
    const Math::MatrixNd & qdot_minus_dirs,
    ConstraintSet   & CS,
    ADConstraintSet * fd_CS,
    Math::VectorNd  & qdot_plus,
    Math::MatrixNd  & fd_qdot_plus
    );

//RigidBodyDynamics::FD::ForwardDynamicsConstraintsDirect(RigidBodyDynamics::Model&,
//     ADModel*,
//     Eigen::Matrix<double, -1, 1, 0, -1, 1> const&,
//     Eigen::Matrix<double, -1, -1, 0, -1, -1> const&,
//     Eigen::Matrix<double, -1, 1, 0, -1, 1> const&,
//     Eigen::Matrix<double, -1, -1, 0, -1, -1> const&,
//     Eigen::Matrix<double, -1, 1, 0, -1, 1> const&,
//     Eigen::Matrix<double, -1, -1, 0, -1, -1> const&,
//     RigidBodyDynamics::ConstraintSet&,
//     RigidBodyDynamics::ADConstraintSet*,
//     Eigen::Matrix<double, -1, 1, 0, -1, 1>&,
//     Eigen::Matrix<double, -1, -1, 0, -1, -1>&)

RBDL_DLLAPI
void ForwardDynamicsConstraintsDirect (
    Model   &model,
    ADModel *fd_model, // NULL means execution without fd_model update
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    const Math::VectorNd &qdot,
    const Math::MatrixNd &qdot_dirs,
    const Math::VectorNd &tau,
    const Math::MatrixNd &tau_dirs,
    ConstraintSet   &cs,
    ADConstraintSet *fd_cs,
    Math::VectorNd  &qddot,
    Math::MatrixNd  &fd_qddot
    );


RBDL_DLLAPI
void ComputeConstraintImpulsesDirect (
    Model & model,
    ADModel * fd_model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot_minus,
    const Math::MatrixNd & qdot_minus_dirs,
    ConstraintSet & CS,
    ADConstraintSet * fd_CS,
    Math::MatrixNd & fd_b,
    std::vector<Math::MatrixNd> & fd_A,
    Math::VectorNd & qdot_plus,
    Math::MatrixNd & fd_qdot_plus
    );


RBDL_DLLAPI
void SolveConstrainedSystemDirect (
    const Math::MatrixNd &H,
    const std::vector<Math::MatrixNd> & H_dirs,
    const Math::MatrixNd &G,
    const std::vector<Math::MatrixNd> & G_dirs,
    const Math::VectorNd & c,
    const Math::MatrixNd & c_dirs,
    const Math::VectorNd & gamma,
    const Math::MatrixNd & gamma_dirs,
    Math::MatrixNd & A,
    std::vector<Math::MatrixNd> & A_dirs,
    Math::VectorNd & b,
    Math::MatrixNd & b_dirs,
    Math::VectorNd & x,
    Math::MatrixNd & x_fd,
    Math::LinearSolver & linear_solver,
    int ndirs
);

/** \brief Computes the effect of external forces on the generalized accelerations.
 *
 * This function is essentially similar to ForwardDynamics() except that it
 * tries to only perform computations of variables that change due to
 * external forces defined in f_t.
 */
RBDL_DLLAPI
void ForwardDynamicsAccelerationDeltas (
    Model &model,
    ADModel *fd_model,
    ConstraintSet &CS,
    ADConstraintSet &ad_CS,
    Math::VectorNd &QDDot_t,
    Math::MatrixNd &ad_QDDot_t,
    const unsigned int body_id,
    const std::vector<Math::SpatialVector> &f_t,
    const std::vector<Math::MatrixNd> &ad_f_t
);

/** \brief Compute only the effects of external forces on the generalized accelerations
 *
 * This function is a reduced version of ForwardDynamics() which only
 * computes the effects of the external forces on the generalized
 * accelerations.
 *
 */
RBDL_DLLAPI
void ForwardDynamicsApplyConstraintForces (
    Model &model,
    ADModel *fd_model,
    const Math::VectorNd &Tau,
    const Math::MatrixNd &Tau_dirs,
    ConstraintSet &CS,
    ADConstraintSet &ad_CS,
    Math::VectorNd &QDDot,
    Math::MatrixNd &ad_QDDot
);

RBDL_DLLAPI
void ForwardDynamicsContactsKokkevis (
  Model &model,
  ADModel *fd_model,
  const Math::VectorNd &q,
  const Math::MatrixNd &q_dirs,
  const Math::VectorNd &qdot,
  const Math::MatrixNd &qdot_dirs,
  const Math::VectorNd &tau,
  const Math::MatrixNd &tau_dirs,
  ConstraintSet &CS,
  ADConstraintSet &fd_cs,
  Math::VectorNd &qddot,
  Math::MatrixNd  &fd_qddot
);

// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

/* RBDL_CONTACTS_FD_H */
#endif
