#ifndef RBDL_CONTACTS_ED_H
#define RBDL_CONTACTS_ED_H

#include <rbdl/rbdl_math.h>
#include <rbdl/rbdl_mathutils.h>

#include "rbdl/Constraints.h"

#include "ModelED.h"
// #include "KinematicsAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

struct EDConstraintSet {
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

  std::vector<Math::SpatialTransform> prev_body_X_1_dirs;
  std::vector<Math::SpatialTransform> prev_body_X_2_dirs;

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

  EDConstraintSet() : CS (0) {}
  EDConstraintSet(const ConstraintSet &CS, int dof_count);

  ~EDConstraintSet() {
    // this class does not own this pointer
    CS = nullptr;
  }

  void resize_directions (const unsigned int& requested_ndirs);
};

// -----------------------------------------------------------------------------
namespace ED {
// -----------------------------------------------------------------------------

RBDL_DLLAPI void CalcConstraintsPositionError (
    Model &model,
    EDModel &ed_model,
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    ConstraintSet &cs,
    EDConstraintSet &ed_cs,
    Math::VectorNd &err,
    Math::MatrixNd &ed_err,
    bool update_kinematics
    );

RBDL_DLLAPI
void CalcConstraintsJacobian (
    Model   & model,
    EDModel & ed_model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    ConstraintSet  & CS,
    EDConstraintSet & ed_CS,
    Math::MatrixNd  & G,
    std::vector<Math::MatrixNd> & G_dirs,
    bool update_kinematics = true
    );

RBDL_DLLAPI
void CalcConstrainedSystemVariables (
    Model &model,
    EDModel &ed_model,
    const Math::VectorNd  &q,
    const Math::MatrixNd  &q_dirs,
    const Math::VectorNd  &qdot,
    const Math::MatrixNd  &qdot_dirs,
    const Math::VectorNd  &tau,
    const Math::MatrixNd  &tau_dirs,
    ConstraintSet   &CS,
    EDConstraintSet &ed_CS);

RBDL_DLLAPI
void ForwardDynamicsConstraintsDirect (
    Model   & model,
    EDModel & ed_model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot,
    const Math::MatrixNd & qdot_dirs,
    const Math::VectorNd & tau,
    const Math::MatrixNd & tau_dirs,
    ConstraintSet   & CS,
    EDConstraintSet & ed_CS,
    Math::VectorNd  & qddot,
    Math::MatrixNd  & ed_qddot
    );

RBDL_DLLAPI
void ComputeConstraintImpulsesDirect (
    Model & model,
    EDModel & ed_model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot_minus,
    const Math::MatrixNd & qdot_minus_dirs,
    ConstraintSet   & CS,
    EDConstraintSet & ed_CS,
    Math::VectorNd  & qdot_plus,
    Math::MatrixNd  & ed_qdot_plus
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
    Math::MatrixNd & x_ad,
    Math::LinearSolver & linear_solver,
    unsigned ndirs
    );

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

/* RBDL_CONTACTS_ED_H */
#endif
