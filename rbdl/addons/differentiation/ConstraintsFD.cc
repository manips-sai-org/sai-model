#include "ConstraintsFD.h"
#include "FdModelEntry.h"

using namespace std;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;

RBDL_DLLAPI void CalcConstrainedSystemVariables (
    Model &model,
    ADModel *fd_model, // NULL means execution without fd_model update
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs,
    const VectorNd &tau,
    const MatrixNd &tau_dirs,
    ConstraintSet   &cs,
    ADConstraintSet *fd_cs) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());

  double const h = 1e-8;

  CalcConstrainedSystemVariables(model, q, qdot, tau, cs);

  ConstraintSet * csh;

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model *modelh = &model;
    VectorNd qh    = q + h * q_dirs.col(idir);
    VectorNd qdoth = qdot + h * qdot_dirs.col(idir);
    VectorNd tauh  = tau + h * tau_dirs.col(idir);

    if (fd_model) {
      modelh = new Model(model);
    }

    if (fd_cs) {
      csh = new ConstraintSet(cs);
    } else {
      csh = &cs;
    }

    CalcConstrainedSystemVariables(*modelh, qh, qdoth, tauh, *csh);

    // TODO make it analogue to ADMOdel with pointer arithmetic
    if(fd_cs) {
      computeFDEntry(cs, *csh, h, idir, *fd_cs);
    }

    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
}


RBDL_DLLAPI void ForwardDynamicsConstraintsDirect (
    Model &model,
    ADModel *fd_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs,
    const VectorNd &tau,
    const MatrixNd &tau_dirs,
    ConstraintSet   &cs,
    ADConstraintSet *fd_cs,
    VectorNd  &qddot,
    MatrixNd  &fd_qddot) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());
  assert(ndirs == fd_qddot.cols());

  double const h = 1e-8;
  ConstraintSet const cs_in = cs;

  ForwardDynamicsConstraintsDirect(model, q, qdot, tau, cs, qddot);

  Model *modelh = &model;
  ConstraintSet csh = cs_in;

  for (unsigned idir = 0; idir < ndirs; idir++) {
    VectorNd qh    = q + h * q_dirs.col(idir);
    VectorNd qdoth = qdot + h * qdot_dirs.col(idir);
    VectorNd tauh  = tau + h * tau_dirs.col(idir);
    VectorNd qddoth = VectorNd::Zero(qddot.rows());

    if (fd_model) {
      modelh = new Model(model);
    }

    ForwardDynamicsConstraintsDirect(*modelh, qh, qdoth, tauh, csh, qddoth);

    fd_qddot.col(idir) = (qddoth - qddot) / h;
    // TODO add pointer arithmetic for fd_cs
    // computeFDEntry(cs, csh, h, idir, fd_cs);

    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
}

RBDL_DLLAPI void CalcConstraintsJacobian(
    Model & model,
    ADModel * fd_model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    ConstraintSet & cs,
    ADConstraintSet & fd_CS,
    MatrixNd & G,
    vector<MatrixNd> & G_dirs) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == G_dirs.size());

  bool const update_kinematics = true;
  double const h = 1e-8;
  ConstraintSet const cs_in = cs;

  CalcConstraintsJacobian(model, q, cs, G, update_kinematics);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    VectorNd qh = q + h * q_dirs.col(idir);
    ConstraintSet csh = cs_in;
    MatrixNd Gh = MatrixNd::Zero (G.rows(), G.cols());

    Model *modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }

    CalcConstraintsJacobian(*modelh, qh, csh, Gh, update_kinematics);
    G_dirs[idir] = (Gh - G) / h;

    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
}

RBDL_DLLAPI void ComputeConstraintImpulsesDirect (
    Model & model,
    ADModel * fd_model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    const VectorNd & qdot_minus,
    const MatrixNd & qdot_minus_dirs,
    ConstraintSet & cs,
    ADConstraintSet * fd_cs,
    VectorNd & qdot_plus,
    MatrixNd & fd_qdot_plus
) {
  unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_minus_dirs.cols());
  assert(ndirs == fd_qdot_plus.cols());

  double h = 1e-8;
  ConstraintSet cs_in = cs;
  ComputeConstraintImpulsesDirect(model, q, qdot_minus, cs, qdot_plus);
  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }
    ConstraintSet csh = cs_in;

    VectorNd qh  = q + h * q_dirs.col(idir);
    VectorNd qdh = qdot_minus + h * qdot_minus_dirs.col(idir);
    VectorNd qdot_plush(model.dof_count);
    ComputeConstraintImpulsesDirect(*modelh, qh, qdh, csh, qdot_plush);
    fd_qdot_plus.col(idir) = (qdot_plush - qdot_plus) / h;

    if (fd_cs) {
      computeFDEntry(cs, csh, h, idir, *fd_cs);
    }
    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
}

RBDL_DLLAPI
void ComputeConstraintImpulsesDirect (
    Model & model,
    ADModel * fd_model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    const VectorNd & qdot_minus,
    const MatrixNd & qdot_minus_dirs,
    ConstraintSet & cs,
    ADConstraintSet * fd_cs,
    MatrixNd & fd_b,
    vector<MatrixNd> & fd_A,
    VectorNd & qdot_plus,
    MatrixNd & fd_qdot_plus
) {
  unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_minus_dirs.cols());
  assert(ndirs == fd_qdot_plus.cols());

  double h = 1e-8;
  ConstraintSet cs_in = cs;
  ComputeConstraintImpulsesDirect(model, q, qdot_minus, cs, qdot_plus);
  VectorNd b_ref = cs.b;
  MatrixNd A_ref = cs.A;
  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }
    ConstraintSet csh = cs_in;

    VectorNd qh  = q + h * q_dirs.col(idir);
    VectorNd qdh = qdot_minus + h * qdot_minus_dirs.col(idir);
    VectorNd qdot_plush(model.dof_count);
    ComputeConstraintImpulsesDirect(*modelh, qh, qdh, csh, qdot_plush);
    fd_qdot_plus.col(idir) = (qdot_plush - qdot_plus) / h;
    fd_b.col(idir) = (cs.b - b_ref) / h;
    fd_A[idir] = (cs.A - A_ref) / h;

    if (fd_cs) {
      computeFDEntry(cs, csh, h, idir, *fd_cs);
    }
    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
}

RBDL_DLLAPI
void SolveConstrainedSystemDirect (
    const MatrixNd &H,
    const vector<MatrixNd> & H_dirs,
    const MatrixNd &G,
    const vector<MatrixNd> & G_dirs,
    const VectorNd & c,
    const MatrixNd & c_dirs,
    const VectorNd & gamma,
    const MatrixNd & gamma_dirs,
    MatrixNd & A,
    vector<MatrixNd> & A_dirs,
    VectorNd & b,
    MatrixNd & b_dirs,
    VectorNd & x,
    MatrixNd & x_fd,
    LinearSolver & linear_solver,
    int ndirs
) {
  double h = 1e-8;

  VectorNd qddot(H.rows());
  VectorNd lambda(H.rows());

  RigidBodyDynamics::SolveConstrainedSystemDirect(
      H, G, c, gamma, qddot, lambda,
      A, b, x, linear_solver);
  for (int i = 0; i < ndirs; i++) {
    MatrixNd Ah(A);
    VectorNd bh(b);
    VectorNd xh(x);
    RigidBodyDynamics::SolveConstrainedSystemDirect (H + h * H_dirs[i],
        G + h * G_dirs[i], c + h * c_dirs.col(i), gamma + h * gamma_dirs.col(i),
        qddot, lambda, Ah, bh, xh, linear_solver);
    A_dirs[i] = (Ah - A) / h;
    b_dirs.col(i) = (bh - b) / h;
    x_fd.col(i) = (xh - x) / h;
  }
}

RBDL_DLLAPI
void ForwardDynamicsAccelerationDeltas (
  Model &model, ADModel *fd_model,
  ConstraintSet &cs, ADConstraintSet &fd_cs,
  VectorNd &QDDot_t, MatrixNd &fd_QDDot_t,
  const unsigned int body_id,
  const std::vector<SpatialVector> &f_t,
  const std::vector<MatrixNd> &f_t_dirs
) {
  const unsigned int ndirs = f_t_dirs[0].cols();

  double const h = 1e-8;
  ConstraintSet const cs_in = cs;

  ForwardDynamicsAccelerationDeltas (model, cs, QDDot_t, body_id, f_t);

  std::vector<SpatialVector> f_t_h (f_t.size(), SpatialVector::Zero());
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < f_t.size(); ++i)
    {
      f_t_h[i] = f_t[i] + h * f_t_dirs[i].col(idir);
    }

    Model *modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }
    ConstraintSet csh = cs_in;

    VectorNd QDDot_t_h (QDDot_t);
    ForwardDynamicsAccelerationDeltas (model, csh, QDDot_t_h, body_id, f_t_h);

    fd_QDDot_t.col(idir) = (QDDot_t_h - QDDot_t) / h;
    // TODO add pointer arithmetic for fd_cs
    computeFDEntry(cs, csh, h, idir, fd_cs);

    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
}

RBDL_DLLAPI
void ForwardDynamicsApplyConstraintForces (
    Model &model,
    ADModel *fd_model,
    const VectorNd &tau,
    const MatrixNd &tau_dirs,
    ConstraintSet &cs,
    ADConstraintSet &ad_cs,
    VectorNd &qddot,
    MatrixNd &fd_qddot
) {
  unsigned const ndirs = tau_dirs.cols();
  assert(ndirs == tau_dirs.cols());
  assert(ndirs == fd_qddot.cols());

  double const h = 1e-8;
  ConstraintSet const cs_in = cs;

  ForwardDynamicsApplyConstraintForces (model, tau, cs, qddot);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    VectorNd tauh  = tau + h * tau_dirs.col(idir);
    VectorNd qddoth (qddot);

    Model *modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }
    ConstraintSet csh = cs_in;

    ForwardDynamicsApplyConstraintForces (*modelh, tauh, cs, qddoth);

    fd_qddot.col(idir) = (qddoth - qddot) / h;
    // TODO add pointer arithmetic for fd_cs
    // computeFDEntry(cs, csh, h, idir, fd_cs);

    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
}

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
  ConstraintSet &cs,
  ADConstraintSet &fd_cs,
  Math::VectorNd &qddot,
  Math::MatrixNd &fd_qddot
) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());
  assert(ndirs == fd_qddot.cols());

  double const h = 1e-8;
  ConstraintSet const cs_in = cs;

  ForwardDynamicsContactsKokkevis (model, q, qdot, tau, cs, qddot);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model *modelh = &model;
    VectorNd qh    = q + h * q_dirs.col(idir);
    VectorNd qdoth = qdot + h * qdot_dirs.col(idir);
    VectorNd tauh  = tau + h * tau_dirs.col(idir);
    ConstraintSet csh = cs_in;
    VectorNd qddoth = VectorNd::Zero(qddot.rows());

    if (fd_model) {
      modelh = new Model(model);
    }

    ForwardDynamicsContactsKokkevis (*modelh, qh, qdoth, tauh, csh, qddoth);

    fd_qddot.col(idir) = (qddoth - qddot) / h;
    // TODO add pointer arithmetic for fd_cs
    // computeFDEntry(cs, csh, h, idir, fd_cs);

    if (fd_model) {
      computeFDEntry(model, *modelh, h, idir, *fd_model);
      delete modelh;
    }
  }
}

// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
