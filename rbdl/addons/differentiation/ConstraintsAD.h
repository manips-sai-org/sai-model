#ifndef RBDL_CONTACTS_AD_H
#define RBDL_CONTACTS_AD_H

#include <rbdl/rbdl_math.h>
#include <rbdl/rbdl_mathutils.h>

#include "rbdl/Constraints.h"

#include "ModelAD.h"
#include "KinematicsAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

struct ADConstraintSet {
  unsigned int ndirs;
  unsigned int dof_count;

  const ConstraintSet* CS;

  std::vector<Math::MatrixNd> G;
  std::vector<Math::MatrixNd> Gi;
  std::vector<Math::MatrixNd> GSpi;
  std::vector<Math::MatrixNd> GSsi;
  std::vector<Math::MatrixNd> GSJ;

  std::vector<Math::MatrixNd> A;
  std::vector<Math::MatrixNd> H;
  std::vector<Math::MatrixNd> K;


  Math::MatrixNd a;
  Math::MatrixNd b;
  Math::MatrixNd v_plus;
  Math::MatrixNd x;
  Math::MatrixNd impulse;
  Math::MatrixNd QDDot_t;
  Math::MatrixNd QDDot_0;
  Math::MatrixNd C;
  Math::MatrixNd gamma;
  Math::MatrixNd force;
  Math::MatrixNd err;
  Math::MatrixNd errd;

  std::vector<Math::MatrixNd> f_t;
  std::vector<Math::MatrixNd> f_ext_constraints;
  std::vector<Math::MatrixNd> point_accel_0;

  std::vector<Math::MatrixNd> d_pA;
  std::vector<Math::MatrixNd> d_a;
  Math::MatrixNd d_u;
  std::vector<Math::MatrixNd> d_multdof3_u;

  Math::MatrixNd point_accel_t;
  Math::MatrixNd point_global;

  ADConstraintSet();
  ADConstraintSet(const ConstraintSet &CS, int dof_count);

  ~ADConstraintSet() {
    // set reference to NULL
    CS = nullptr;
  }

  void resize_directions (const unsigned int& requested_ndirs);
};

// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI void CalcConstraintsPositionError (
    Model &model,
    ADModel &ad_model,
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    ConstraintSet &cs,
    ADConstraintSet &ad_cs,
    Math::VectorNd &err,
    Math::MatrixNd &ad_err,
    bool update_kinematics
    );

RBDL_DLLAPI
void CalcConstraintsJacobian(
		Model   & model,
		ADModel & ad_model,
		const Math::VectorNd & q,
		const Math::MatrixNd & q_dirs,
    ConstraintSet  & CS,
		ADConstraintSet & ad_CS,
		Math::MatrixNd  & G,
		std::vector<Math::MatrixNd> & G_dirs,
		bool update_kinematics = true
		);

RBDL_DLLAPI void CalcConstrainedSystemVariables (
    Model &model,
    ADModel &ad_model,
    const Math::VectorNd  &q,
    const Math::MatrixNd  &q_dirs,
    const Math::VectorNd  &qdot,
    const Math::MatrixNd  &qdot_dirs,
    const Math::VectorNd  &tau,
    const Math::MatrixNd  &tau_dirs,
    ConstraintSet   &CS,
    ADConstraintSet &ad_CS);

RBDL_DLLAPI
void ForwardDynamicsConstraintsDirect (
    Model   & model,
    ADModel & ad_model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot,
    const Math::MatrixNd & qdot_dirs,
    const Math::VectorNd & tau,
    const Math::MatrixNd & tau_dirs,
    ConstraintSet   & CS,
    ADConstraintSet & ad_CS,
    Math::VectorNd  & qddot,
    Math::MatrixNd  & ad_qddot
    );

/*
RBDL_DLLAPI
void ForwardDynamicsContactsRangeSpaceSparse (
				Model &model,
				const Math::VectorNd &Q,
				const Math::VectorNd &QDot,
				const Math::VectorNd &Tau,
				ConstraintSet &CS,
				Math::VectorNd &QDDot
				);

RBDL_DLLAPI
void ForwardDynamicsContactsNullSpace (
				Model &model,
				const Math::VectorNd &Q,
				const Math::VectorNd &QDot,
				const Math::VectorNd &Tau,
				ConstraintSet &CS,
				Math::VectorNd &QDDot
				);

*/

RBDL_DLLAPI
void ComputeConstraintImpulsesDirect(
		Model & model,
		ADModel & ad_model,
		const Math::VectorNd & q,
		const Math::MatrixNd & q_dirs,
		const Math::VectorNd & qdot_minus,
		const Math::MatrixNd & qdot_minus_dirs,
		ConstraintSet   & CS,
		ADConstraintSet & ad_CS,
		Math::VectorNd  & qdot_plus,
		Math::MatrixNd  & ad_qdot_plus
		);

/*
RBDL_DLLAPI
void ComputeContactImpulsesRangeSpaceSparse (
				Model &model,
				const Math::VectorNd &Q,
				const Math::VectorNd &QDotMinus,
				ConstraintSet &CS,
				Math::VectorNd &QDotPlus
				);

RBDL_DLLAPI
void ComputeContactImpulsesNullSpace (
				Model &model,
				const Math::VectorNd &Q,
				const Math::VectorNd &QDotMinus,
				ConstraintSet &CS,
				Math::VectorNd &QDotPlus
				);
*/

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
    Math::MatrixNd & x_ad,
    Math::LinearSolver & linear_solver,
    unsigned ndirs
    );

/*
RBDL_DLLAPI
void SolveContactSystemRangeSpaceSparse (
				Model &model,
				Math::MatrixNd &H,
				const Math::MatrixNd &G,
				const Math::VectorNd &c,
				const Math::VectorNd &gamma,
				Math::VectorNd &qddot,
				Math::VectorNd &lambda,
				Math::MatrixNd &K,
				Math::VectorNd &a,
				Math::LinearSolver linear_solver
				);

RBDL_DLLAPI
void SolveContactSystemNullSpace (
				Math::MatrixNd &H,
				const Math::MatrixNd &G,
				const Math::VectorNd &c,
				const Math::VectorNd &gamma,
				Math::VectorNd &qddot,
				Math::VectorNd &lambda,
				Math::MatrixNd &Y,
				Math::MatrixNd &Z,
				Math::VectorNd &qddot_y,
				Math::VectorNd &qddot_z,
				Math::LinearSolver &linear_solver
				);
*/

/** \brief Computes the effect of external forces on the generalized accelerations.
 *
 * This function is essentially similar to ForwardDynamics() except that it
 * tries to only perform computations of variables that change due to
 * external forces defined in f_t.
 */
RBDL_DLLAPI
void ForwardDynamicsAccelerationDeltas (
    Model &model,
    ADModel &ad_model,
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
    ADModel &ad_model,
    const Math::VectorNd &Tau,
    const Math::MatrixNd &Tau_dirs,
    ConstraintSet &CS,
    ADConstraintSet &ad_CS,
    Math::VectorNd &QDDot,
    Math::MatrixNd &ad_QDDot
);

RBDL_DLLAPI
void ForwardDynamicsContactsKokkevis (
  Model   & model,
  ADModel & ad_model,
  const Math::VectorNd & q,
  const Math::MatrixNd & q_dirs,
  const Math::VectorNd & qdot,
  const Math::MatrixNd & qdot_dirs,
  const Math::VectorNd & tau,
  const Math::MatrixNd & tau_dirs,
  ConstraintSet   & CS,
  ADConstraintSet & ad_CS,
  Math::VectorNd  & qddot,
  Math::MatrixNd  & ad_qddot
);

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

/* RBDL_CONTACTS_AD_H */
#endif
