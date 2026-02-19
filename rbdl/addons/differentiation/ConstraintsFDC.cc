#include "ConstraintsFDC.h"
#include "ModelEntryFDC.h"

#include <cfloat>
#include <cmath>

using namespace std;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FDC {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;

// define constant difference perturbation according to theory,
// i.e. EPS = eps^(1/3) = cubic_root(eps).
// Here, we us cmath cubic root implementation, i.e. cubic_root = cbrt.
// NOTE assumes that directions are normalized to one
const double EPS = cbrt(DBL_EPSILON);
const double EPSx2 = 2*EPS; // for convenience, we directly compute the denominator

RBDL_DLLAPI
void CalcConstrainedSystemVariables (
    Model &model,
    ADModel *fd_model, // NULL means execution without fd_model update
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs,
    const VectorNd &tau,
    const MatrixNd &tau_dirs,
    ConstraintSet  &cs,
    ADConstraintSet *fd_cs // NULL means execution without fd_model update
) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }

    ConstraintSet * csh;
    if (fd_cs) {
      csh = new ConstraintSet(cs);
    } else {
      csh = &cs;
    }

    // forward perturbation
    CalcConstrainedSystemVariables(
      *modelh,
      q + EPS * q_dirs.col(idir),
      qdot + EPS * qdot_dirs.col(idir),
      tau + EPS * tau_dirs.col(idir),
      *csh
    );

    // backward perturbation
    CalcConstrainedSystemVariables(
      model,
      q - EPS * q_dirs.col(idir),
      qdot - EPS * qdot_dirs.col(idir),
      tau - EPS * tau_dirs.col(idir),
      cs
    );

    if (fd_model) {
      computeFDCEntry(*modelh, model, EPS, idir, *fd_model);
      delete modelh;
    }

    if (fd_cs) {
      computeFDCEntry(*csh, cs, EPS, idir, *fd_cs);
      delete csh;
    }
  }

  // nominal evaluation
  CalcConstrainedSystemVariables(model, q, qdot, tau, cs);

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
  MatrixNd  &fd_qddot
) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());
  assert(ndirs == fd_qddot.cols());

  ConstraintSet const cs_in = cs;


  // temporary variables
  VectorNd qddotph (qddot);
  VectorNd qddotmh (qddot);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }

    ConstraintSet * csh;
    if (fd_cs) {
      csh = new ConstraintSet(cs);
    } else {
      csh = &cs;
    }

    ForwardDynamicsConstraintsDirect(
      *modelh,
      q + EPS * q_dirs.col(idir),
      qdot + EPS * qdot_dirs.col(idir),
      tau + EPS * tau_dirs.col(idir),
      *csh,
      qddotph
    );

    ForwardDynamicsConstraintsDirect(
      model,
      q - EPS * q_dirs.col(idir),
      qdot - EPS * qdot_dirs.col(idir),
      tau - EPS * tau_dirs.col(idir),
      cs,
      qddotmh
    );

    fd_qddot.col(idir) = (qddotph - qddotmh) / EPSx2;

    if (fd_model) {
      computeFDCEntry(*modelh, model, EPS, idir, *fd_model);
      delete modelh;
    }

    if (fd_cs) {
      computeFDCEntry(*csh, cs, EPS, idir, *fd_cs);
      delete csh;
    }
  }

  ForwardDynamicsConstraintsDirect(model, q, qdot, tau, cs, qddot);
}

RBDL_DLLAPI void ForwardDynamicsContactsKokkevis (
  Model &model,
  ADModel *fd_model,
  const VectorNd &q,
  const MatrixNd &q_dirs,
  const VectorNd &qdot,
  const MatrixNd &qdot_dirs,
  const VectorNd &tau,
  const MatrixNd &tau_dirs,
  ConstraintSet   &cs,
  ADConstraintSet &fd_cs,
  VectorNd  &qddot,
  MatrixNd  &fd_qddot
) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());
  assert(ndirs == fd_qddot.cols());

  ConstraintSet const cs_in = cs;
  ForwardDynamicsContactsKokkevis(model, q, qdot, tau, cs, qddot);

  // temporary variables
  VectorNd qddotph (qddot);
  VectorNd qddotmh (qddot);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelhp;
    Model * modelhm;
    ConstraintSet * cshp;
    ConstraintSet * cshm;
    if (fd_model) {

    } else {
      modelhp = &model;
      modelhm = &model;
      cshp = &cs;
      cshm = &cs;
    }
    cshp = new ConstraintSet(cs_in);
    cshm = new ConstraintSet(cs_in);
    modelhp = new Model(model);
    modelhm = new Model(model);


    ForwardDynamicsContactsKokkevis(
      *modelhp,
      q + EPS * q_dirs.col(idir),
      qdot + EPS * qdot_dirs.col(idir),
      tau + EPS * tau_dirs.col(idir),
      *cshp,
      qddotph
    );

    ForwardDynamicsContactsKokkevis(
      *modelhm,
      q - EPS * q_dirs.col(idir),
      qdot - EPS * qdot_dirs.col(idir),
      tau - EPS * tau_dirs.col(idir),
      *cshm,
      qddotmh
    );

    fd_qddot.col(idir) = (qddotph - qddotmh) / EPSx2;

    // TODO add pointer arithmetic for fd_cs
    // computeFDEntry(cs, csh, EPS, idir, fd_cs);

    if (fd_model) {
      computeFDCEntry(*modelhp, *modelhm, EPS, idir, *fd_model);
      computeFDCEntry(*cshp, *cshm, EPS, idir, fd_cs);
    }
    delete cshp;
    delete cshm;

    delete modelhp;
    delete modelhm;
  }
}

RBDL_DLLAPI void CalcConstraintsJacobian(
    Model & model,
    ADModel * fd_model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    ConstraintSet & cs,
    ADConstraintSet * fd_cs,
    MatrixNd & G,
    vector<MatrixNd> & G_dirs
) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == G_dirs.size());

  bool const update_kinematics = true;

  for (unsigned idir = 0; idir < ndirs; idir++) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }

    ConstraintSet * csh;
    if (fd_cs) {
      csh = new ConstraintSet(cs);
    } else {
      csh = &cs;
    }

    MatrixNd Ghp = MatrixNd::Zero (G.rows(), G.cols());
    MatrixNd Ghm = MatrixNd::Zero (G.rows(), G.cols());

    CalcConstraintsJacobian(
      *modelh,
      q + EPS * q_dirs.col(idir),
      *csh,
      Ghp,
      update_kinematics
    );

    CalcConstraintsJacobian(
      model,
      q - EPS * q_dirs.col(idir),
      cs,
      Ghm,
      update_kinematics
    );

    G_dirs[idir] = (Ghp - Ghm) / EPSx2;

    if (fd_model) {
      computeFDCEntry(*modelh, model, EPS, idir, *fd_model);
      delete modelh;
    }

    if (fd_cs) {
      computeFDCEntry(*csh, cs, EPS, idir, *fd_cs);
      delete csh;
    }
  }

  CalcConstraintsJacobian(model, q, cs, G, update_kinematics);
}

/*
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

    VectorNd qh  = q + EPS * q_dirs.col(idir);
    VectorNd qdh = qdot_minus + EPS * qdot_minus_dirs.col(idir);
    VectorNd qdot_plush(model.dof_count);
    ComputeConstraintImpulsesDirect(*modelh, qh, qdh, csh, qdot_plush);
    fd_qdot_plus.col(idir) = (qdot_plush - qdot_plus) / EPS;

    if (fd_cs) {
      computeFDCEntry(cs, csh, EPS, idir, *fd_cs);
    }
    if (fd_model) {
      computeFDCEntry(model, *modelh, EPS, idir, *fd_model);
      delete modelh;
    }
  }
}
*/
/*
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

    VectorNd qh  = q + EPS * q_dirs.col(idir);
    VectorNd qdh = qdot_minus + EPS * qdot_minus_dirs.col(idir);
    VectorNd qdot_plush(model.dof_count);
    ComputeConstraintImpulsesDirect(*modelh, qh, qdh, csh, qdot_plush);
    fd_qdot_plus.col(idir) = (qdot_plush - qdot_plus) / EPS;
    fd_b.col(idir) = (cs.b - b_ref) / EPS;
    fd_A[idir] = (cs.A - A_ref) / EPS;

    if (fd_cs) {
      computeFDCEntry(cs, csh, EPS, idir, *fd_cs);
    }
    if (fd_model) {
      computeFDCEntry(model, *modelh, EPS, idir, *fd_model);
      delete modelh;
    }
  }
}
*/

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
  VectorNd qddot(H.rows());
  VectorNd lambda(H.rows());

  MatrixNd Ah(A);
  VectorNd bh(b);
  VectorNd xh(x);

  for (int i = 0; i < ndirs; i++) {
    RigidBodyDynamics::SolveConstrainedSystemDirect (
        H + EPS * H_dirs[i],
        G + EPS * G_dirs[i],
        c + EPS * c_dirs.col(i),
        gamma + EPS * gamma_dirs.col(i),
        qddot, lambda,
        Ah, bh, xh,
        linear_solver
    );

    RigidBodyDynamics::SolveConstrainedSystemDirect (
        H - EPS * H_dirs[i],
        G - EPS * G_dirs[i],
        c - EPS * c_dirs.col(i),
        gamma - EPS * gamma_dirs.col(i),
        qddot, lambda,
        A, b, x,
        linear_solver
    );

    A_dirs[i] = (Ah - A) / EPSx2;
    b_dirs.col(i) = (bh - b) / EPSx2;
    x_fd.col(i) = (xh - x) / EPSx2;
  }

  RigidBodyDynamics::SolveConstrainedSystemDirect(
      H, G, c, gamma, qddot, lambda,
      A, b, x, linear_solver);

}

// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
