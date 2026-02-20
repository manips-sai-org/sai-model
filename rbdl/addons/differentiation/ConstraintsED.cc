#include "rbdl/Constraints.h"
#include "ConstraintsED.h"
#include "KinematicsED.h"
#include "DynamicsED.h"
#include "rbdl/SpatialAlgebraOperators.h"
#include "SpatialAlgebraOperatorsAD.h"
#include "rbdl_mathutilsAD.h"

#include <iomanip>

using namespace std;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;

EDConstraintSet::EDConstraintSet(const ConstraintSet & CS_, int _dof_count) :
  CS (0)
{
  // std::cout << "in " << __func__ << std::endl;
  dof_count = _dof_count;
  ndirs = 4 * dof_count;
  // const unsigned int ndirs_ = 4 * dof_count;

  CS = &CS_;

  G.resize(ndirs,  MatrixNd::Zero(CS->G.rows(), CS->G.cols()));
  Gi.resize(ndirs, MatrixNd::Zero(3, dof_count));
  GSpi.resize(ndirs, MatrixNd::Zero(6, dof_count));
  GSsi.resize(ndirs, MatrixNd::Zero(6, dof_count));

  A.resize(ndirs,  MatrixNd::Zero(CS->A.rows(), CS->A.cols()));
  H.resize(ndirs,  MatrixNd::Zero(CS->H.rows(), CS->H.cols()));
  K.resize(ndirs,  MatrixNd::Zero(CS->K.rows(), CS->K.cols()));

  prev_body_X_1_dirs.resize(ndirs, SpatialTransform::Zero());
  prev_body_X_2_dirs.resize(ndirs, SpatialTransform::Zero());

  a             = MatrixNd::Zero(CS->a.rows(), ndirs);
  b             = MatrixNd::Zero(CS->b.rows(), ndirs);
  v_plus        = MatrixNd::Zero(CS->v_plus.rows(), ndirs);
  x             = MatrixNd::Zero(CS->x.rows(), ndirs);
  impulse       = MatrixNd::Zero(CS->impulse.rows(), ndirs);
  QDDot_t       = MatrixNd::Zero(CS->QDDot_t.rows(), ndirs);
  QDDot_0       = MatrixNd::Zero(CS->QDDot_0.rows(), ndirs);
  C             = MatrixNd::Zero(CS->C.rows(), ndirs);
  gamma         = MatrixNd::Zero(CS->gamma.rows(), ndirs);
  force         = MatrixNd::Zero(CS->force.rows(), ndirs);
  err           = MatrixNd::Zero(CS->err.rows(), ndirs);
  errd          = MatrixNd::Zero(CS->errd.rows(), ndirs);

  // Kokkevis values
  f_t.resize(CS->f_t.size(),  MatrixNd::Zero(6, ndirs));
  f_ext_constraints.resize(CS->f_ext_constraints.size(),  MatrixNd::Zero(6, ndirs));
  point_accel_0.resize(CS->point_accel_0.size(), MatrixNd::Zero(3, ndirs));

  d_pA.resize(CS->d_pA.size(),  MatrixNd::Zero(6, ndirs));
  d_a.resize(CS->d_a.size(),  MatrixNd::Zero(6, ndirs));
  d_u = MatrixNd::Zero(CS->d_u.rows(), ndirs);
  d_multdof3_u.resize(CS->d_multdof3_u.size(),  MatrixNd::Zero(3, ndirs));

  point_accel_t = MatrixNd::Zero(3, ndirs);
  point_global = MatrixNd::Zero(3, ndirs);
}

void EDConstraintSet::resize_directions(const unsigned int& requested_ndirs) {
  // std::cout << "in " << __func__ << std::endl;
  QDDot_t       = MatrixNd::Zero(CS->QDDot_t.rows(), requested_ndirs);
  QDDot_0       = MatrixNd::Zero(CS->QDDot_0.rows(), requested_ndirs);

  if (ndirs < requested_ndirs) {
    // std::cout << "resized!" << std::endl;
    ndirs = requested_ndirs;

    G.resize(ndirs,  MatrixNd::Zero(CS->G.rows(), CS->G.cols()));
    Gi.resize(ndirs, MatrixNd::Zero(3, dof_count));
    GSpi.resize(ndirs, MatrixNd::Zero(6, dof_count));
    GSsi.resize(ndirs, MatrixNd::Zero(6, dof_count));

    A.resize(ndirs,  MatrixNd::Zero(CS->A.rows(), CS->A.cols()));
    H.resize(ndirs,  MatrixNd::Zero(CS->H.rows(), CS->H.cols()));
    K.resize(ndirs,  MatrixNd::Zero(CS->K.rows(), CS->K.cols()));

    prev_body_X_1_dirs.resize(ndirs, SpatialTransform::Zero());
    prev_body_X_2_dirs.resize(ndirs, SpatialTransform::Zero());

    a             = MatrixNd::Zero(CS->a.rows(), ndirs);
    b             = MatrixNd::Zero(CS->b.rows(), ndirs);
    v_plus        = MatrixNd::Zero(CS->v_plus.rows(), ndirs);
    x             = MatrixNd::Zero(CS->x.rows(), ndirs);
    impulse       = MatrixNd::Zero(CS->impulse.rows(), ndirs);
    C             = MatrixNd::Zero(CS->C.rows(), ndirs);
    gamma         = MatrixNd::Zero(CS->gamma.rows(), ndirs);
    force         = MatrixNd::Zero(CS->force.rows(), ndirs);
    err           = MatrixNd::Zero(CS->err.rows(), ndirs);
    errd          = MatrixNd::Zero(CS->errd.rows(), ndirs);

    // Kokkevis values
    f_t.resize(CS->f_t.size(),  MatrixNd::Zero(6, ndirs));
    f_ext_constraints.resize(CS->f_t.size(),  MatrixNd::Zero(6, ndirs));
    point_accel_0.resize(CS->point_accel_0.size(), MatrixNd::Zero(3, ndirs));

    d_pA.resize(CS->d_pA.size(),  MatrixNd::Zero(6, ndirs));
    d_a.resize(CS->d_a.size(),  MatrixNd::Zero(6, ndirs));
    d_u = MatrixNd::Zero(3, ndirs);
    d_multdof3_u.resize(CS->d_multdof3_u.size(),  MatrixNd::Zero(3, ndirs));

    point_accel_t = MatrixNd::Zero(3, ndirs);
    point_global = MatrixNd::Zero(3, ndirs);
  }
}

// -----------------------------------------------------------------------------
namespace ED {
// -----------------------------------------------------------------------------

/*
RBDL_DLLAPI
void ComputeConstraintImpulsesDirect (
    Model &model,
    EDModel &ed_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot_minus,
    const MatrixNd &qdot_minus_dirs,
    ConstraintSet &CS,
    EDConstraintSet &ed_CS,
    VectorNd &qdot_plus,
    MatrixNd &ed_qdot_plus) {
  int ndirs = q_dirs.cols();
  assert(ndirs == qdot_minus_dirs.cols());
  assert(ndirs == ed_qdot_plus.cols());

  // Compute H
  UpdateKinematicsCustom(model, ed_model, &q, &q_dirs, 0, 0, 0, 0);
  CompositeRigidBodyAlgorithm(model, ed_model, q, q_dirs, CS.H, ed_CS.H, false);

  // Compute G
  CalcConstraintsJacobian (model, ed_model, q, q_dirs,
                       CS, ed_CS, CS.G, ed_CS.G, false);

  VectorNd c = CS.H * qdot_minus;
  MatrixNd c_ad(CS.H.rows(), ndirs);
  for (int i = 0; i < ndirs; i++) {
    c_ad.block(0, i, CS.H.rows(), 1) = CS.H * qdot_minus_dirs.col(i) + ed_CS.H[i] * qdot_minus;
  }

  SolveConstrainedSystemDirect (CS.H, ed_CS.H, CS.G, ed_CS.G, c, c_ad, CS.v_plus,
                            ed_CS.v_plus, CS.A, ed_CS.A, CS.b, ed_CS.b,
                            CS.x, ed_CS.x, CS.linear_solver, ndirs);

  // Copy back QDotPlus
  // derivative evaluation
  ed_qdot_plus = ed_CS.x.block(0, 0, model.dof_count, ndirs);
  // nominal evaluation
  qdot_plus = CS.x.segment(0, model.dof_count);

  // Copy back constraint impulses
  ed_CS.impulse = ed_CS.x.block(model.dof_count, 0, CS.size(), ndirs);
  CS.impulse = CS.x.segment(model.dof_count, CS.size());
}
*/

RBDL_DLLAPI void SolveConstrainedSystemDirect (
    const MatrixNd &H,
    const vector<MatrixNd> &H_dirs,
    const MatrixNd &G,
    const vector<MatrixNd> &G_dirs,
    const VectorNd &c,
    const MatrixNd &c_dirs,
    const VectorNd &gamma,
    const MatrixNd &gamma_dirs,
    MatrixNd &A,
    vector<MatrixNd> &A_dirs,
    VectorNd &b,
    MatrixNd &b_dirs,
    VectorNd &x,
    MatrixNd &x_ad,
    LinearSolver &linear_solver,
    unsigned ndirs) {
  assert(ndirs <= H_dirs.size());
  assert(ndirs <= G_dirs.size());
  assert(ndirs <= c_dirs.cols());
  assert(ndirs <= gamma_dirs.cols());
  assert(ndirs <= A_dirs.size());
  assert(ndirs <= b_dirs.cols());
  assert(ndirs <= x_ad.cols());

  // derivative construction
  for (unsigned idir = 0; idir < ndirs; idir++) {
    A_dirs[idir].block(0, 0, c.rows(), c.rows())            = H_dirs[idir];
    A_dirs[idir].block(0, c.rows(), c.rows(), gamma.rows()) =
        G_dirs[idir].transpose();
    A_dirs[idir].block(c.rows(), 0, gamma.rows(), c.rows()) = G_dirs[idir];
    b_dirs.block(0, idir, c.rows(), 1)                      = c_dirs.col(idir);
    b_dirs.block(c.rows(), idir, gamma.rows(), 1)           = gamma_dirs.col(idir);
  }
  // nominal construction
  //    Build the system: Copy H
  A.block(0, 0, c.rows(), c.rows()) = H;
  //    Copy G and G^T
  A.block(0, c.rows(), c.rows(), gamma.rows()) = G.transpose();
  A.block(c.rows(), 0, gamma.rows(), c.rows()) = G;
  //    Build the system: Copy -C + \tau
  b.block(0, 0, c.rows(), 1) = c;
  b.block(c.rows(), 0, gamma.rows(), 1) = gamma;

  LOG << "A = " << std::endl << A << std::endl;
  LOG << "b = " << std::endl << b << std::endl;

  switch (linear_solver) {
    case (LinearSolverPartialPivLU) :
#ifdef RBDL_USE_SIMPLE_MATH
      {
        cerr << __FILE__ << " " << __LINE__
             << ": Simple Math not supported." << endl;
      }
#else
      {
        Eigen::PartialPivLU<MatrixNd::PlainObject> A_LU = A.partialPivLu();
        // nominal code
        x = A_LU.solve(b);
        // derivative code
        for (unsigned idir = 0; idir < ndirs; idir++) {
          x_ad.col(idir) = A_LU.solve(b_dirs.col(idir) - A_dirs[idir] * x);
        }
      }
#endif
      break;
    case (LinearSolverColPivHouseholderQR) :
      {
        Eigen::ColPivHouseholderQR<MatrixNd::PlainObject> A_CPQR =
            A.colPivHouseholderQr();
        // nominal code
        x = A_CPQR.solve(b);
        // derivative code
        for (unsigned idir = 0; idir < ndirs; idir++) {
          x_ad.col(idir) = A_CPQR.solve(b_dirs.col(idir) - A_dirs[idir] * x);
        }
      }
      break;
    case (LinearSolverHouseholderQR) :
      {
        Eigen::HouseholderQR<MatrixNd::PlainObject> A_QR = A.householderQr();
        // nominal evaluation
        x = A_QR.solve(b);
        // derivative evaluation
        for (unsigned idir = 0; idir < ndirs; idir++) {
          x_ad.col(idir) = A_QR.solve(b_dirs.col(idir) - A_dirs[idir] * x);
        }
      }
      break;
//    case (LinearSolverLLT) :
//      {
//        Eigen::LLT<MatrixNd::PlainObject> A_LLT = A.llt();
//        // nominal code
//        x = A_LLT.solve(b);
//        // derivative code
//        for (int i = 0; i < ndirs; i++) {
//          x_ad.col(i) = A_LLT.solve(b_dirs.col(i) - A_dirs[i] * x);
//        }
//      }
//      break;
/// TODO: Add support for this solver also.
    default:
      LOG << "Error: Invalid linear solver: " << linear_solver << std::endl;
      cerr << __FILE__ << " " << __LINE__
           << ": Error: Invalid linear solver: " << linear_solver << std::endl;
      assert (0);
      break;
  }
}

/*
RBDL_DLLAPI void CalcConstraintsPositionError (
    Model &model,
    EDModel &ed_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    ConstraintSet &cs,
    EDConstraintSet &ed_cs,
    VectorNd &err,
    MatrixNd &ed_err,
    bool update_kinematics) {
  unsigned const ndirs = q_dirs.cols();
  assert(static_cast<unsigned>(err.size()) == cs.size());
  assert(ndirs <= ed_err.cols());

  if(update_kinematics) {
    AD::UpdateKinematicsCustom(model, ed_model,
                           &q, &q_dirs,
                           NULL, NULL,
                           NULL, NULL);
  }

  for (unsigned int i = 0; i < cs.contactConstraintIndices.size(); i++) {
    unsigned const c = cs.contactConstraintIndices[i];
    ed_err.row(c).setZero();
    err[c] = 0.;
  }

  for (unsigned int i = 0; i < cs.loopConstraintIndices.size(); i++) {
    unsigned const lci = cs.loopConstraintIndices[i];

    MatrixNd ed_pos_p(3, ndirs);
    MatrixNd ed_pos_s(3, ndirs);
    vector<Matrix3d> ed_rot_p(ndirs);
    vector<Matrix3d> ed_rot_s(ndirs);
    vector<Matrix3d> ed_rot_ps(ndirs);
    vector<SpatialVector> ed_d(ndirs);

    // Variables used for computations.
    Vector3d pos_p;
    Vector3d pos_s;
    Matrix3d rot_p;
    Matrix3d rot_s;
    Matrix3d rot_ps;
    SpatialVector d;

    // Constraints computed in the predecessor body frame.

    // Compute the orientation of the two constraint frames.
    rot_p = CalcBodyWorldOrientation (model, ed_model, q, q_dirs,
                                      cs.body_p[lci], ed_rot_p, false);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ed_rot_p[idir] = ed_rot_p[idir].transpose() * cs.X_p[lci].E;
    }
    rot_p = rot_p.transpose() * cs.X_p[lci].E;

    rot_s = CalcBodyWorldOrientation (model, ed_model, q, q_dirs,
                                      cs.body_s[lci], ed_rot_s, false);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ed_rot_s[idir] = ed_rot_s[idir].transpose() * cs.X_p[lci].E;
    }
    rot_s = rot_s.transpose() * cs.X_s[lci].E;

    // Compute the orientation from the predecessor to the successor frame.
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ed_rot_ps[idir] =
          ed_rot_p[idir].transpose() * rot_s
          + rot_p.transpose() * ed_rot_s[idir];
    }
    rot_ps = rot_p.transpose() * rot_s;

    // Compute the position of the two contact points.
    pos_p = CalcBodyToBaseCoordinates(model, ed_model, q, q_dirs,
                                      cs.body_p[lci], cs.X_p[lci].r,
                                      ed_pos_p, false);
    pos_s = CalcBodyToBaseCoordinates(model, ed_model, q, q_dirs,
                                      cs.body_s[lci], cs.X_s[lci].r,
                                      ed_pos_s, false);

    // The first three elemenets represent the rotation error.
    // This formulation is equivalent to u * sin(theta), where u and theta are
    // the angle-axis of rotation from the predecessor to the successor frame.
    // These quantities are expressed in the predecessor frame.
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ed_d[idir](0) = -0.5 * (ed_rot_ps[idir](1,2) - ed_rot_ps[idir](2,1));
      ed_d[idir](1) = -0.5 * (ed_rot_ps[idir](2,0) - ed_rot_ps[idir](0,2));
      ed_d[idir](2) = -0.5 * (ed_rot_ps[idir](0,1) - ed_rot_ps[idir](1,0));
    }
    d[0] = -0.5 * (rot_ps(1,2) - rot_ps(2,1));
    d[1] = -0.5 * (rot_ps(2,0) - rot_ps(0,2));
    d[2] = -0.5 * (rot_ps(0,1) - rot_ps(1,0));

    // The last three elements represent the position error.
    // It is equivalent to the difference in the position of the two
    // constraint points.
    // The distance is projected on the predecessor frame to be consistent
    // with the rotation.
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ed_d[idir].segment<3>(3) =
          ed_rot_p[idir].transpose() * (pos_s - pos_p)
          + rot_p.transpose() * (ed_pos_s.col(idir) - ed_pos_p.col(idir));
    }
    d.segment<3>(3) = rot_p.transpose() * (pos_s - pos_p);

    // Project the error on the constraint axis to find the actual error.
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ed_err.col(idir) = cs.constraintAxis[lci].transpose() * ed_d[idir];
    }
    err[lci] = cs.constraintAxis[lci].transpose() * d;
  }
}
*/

RBDL_DLLAPI void CalcConstraintsJacobian(
    Model &model,
    EDModel &ed_model,
    const VectorNd   &q,
    const MatrixNd   &q_dirs,
    ConstraintSet    &CS,
    EDConstraintSet  &ed_CS,
    MatrixNd         &G,
    vector<MatrixNd> &G_dirs,
    bool update_kinematics
) {
  const unsigned int ndirs = q_dirs.cols();
  assert(ndirs <= G_dirs.size());

  // resize if required
  ed_model.resize_directions(ndirs);
  ed_CS.resize_directions(ndirs);

  if (update_kinematics) {
    UpdateKinematicsCustom (
      model, ed_model,
      &q, &q_dirs,
      nullptr, nullptr,
      nullptr, nullptr
    );
  }

  // variables to check whether we need to recompute G.
  ConstraintSet::ConstraintType prev_constraint_type
    = ConstraintSet::ConstraintTypeLast;
  unsigned int prev_body_id_1 = 0;
  // unsigned int prev_body_id_2 = 0;
  SpatialTransform prev_body_X_1;
  SpatialTransform prev_body_X_2;

  for (unsigned int i = 0; i < CS.contactConstraintIndices.size(); i++)
  {
    const unsigned int c = CS.contactConstraintIndices[i];

    // only compute the matrix Gi if actually needed
    if (prev_constraint_type != CS.constraintType[c]
        || prev_body_id_1 != CS.body[c]
        || prev_body_X_1.r != CS.point[c])
    {

      // Compute the jacobian for the point.
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ed_CS.Gi[idir].setZero();
      }
      CS.Gi.setZero();
      CalcPointJacobian (
        model, ed_model, q, q_dirs,
        CS.body[c], CS.point[c],
        CS.Gi, ed_CS.Gi,
        false
      );
      prev_constraint_type = ConstraintSet::ContactConstraint;

      // Update variables for optimization check.
      prev_body_id_1 = CS.body[c];
      // derivative evaluation
      // for (unsigned idir = 0; idir < ndirs; idir++) {
      //   ed_model.prev_body_X_1_dirs[idir].setZero();
      // }
      // nominal evaluation
      prev_body_X_1 = Xtrans(CS.point[c]);
    }

    for(unsigned int j = 0; j < model.dof_count; j++) {
      // derivative evaluation
      MatrixNd gaxis_dirs(3, ndirs);
      for (unsigned idir = 0; idir < ndirs; idir++) {
        gaxis_dirs.col(idir) = ed_CS.Gi[idir].col(j);
        G_dirs[idir](c,j) = gaxis_dirs.col(idir).transpose() * CS.normal[c];
      }
      // nominal evaluation
      Vector3d gaxis (CS.Gi(0,j), CS.Gi(1,j), CS.Gi(2,j));
      G(c,j) = gaxis.transpose() * CS.normal[c];
    }

  }

  // Variables used for computations.
  Vector3d normal;
  SpatialVector axis;
  Vector3d pos_p;
  Matrix3d rot_p;
  SpatialTransform X_0p;

  for (unsigned int i = 0; i < CS.loopConstraintIndices.size(); i++) {
      cerr << __FILE__ << " " << __LINE__
           << ": Error: Loop constraints not supported!" << std::endl;
      assert (0);

  /*
    const unsigned int c = CS.loopConstraintIndices[i];
    // Only recompute variables if necessary.
    if( prev_body_id_1 != CS.body_p[c]
        || prev_body_id_2 != CS.body_s[c]
        || prev_body_X_1.r != CS.X_p[c].r
        || prev_body_X_2.r != CS.X_s[c].r
        || prev_body_X_1.E != CS.X_p[c].E
        || prev_body_X_2.E != CS.X_s[c].E) {

      // Compute the 6D jacobians of the two contact points.
      CS.GSpi.setZero();
      CS.GSsi.setZero();
      CalcPointJacobian6D(model, Q, CS.body_p[c], CS.X_p[c].r, CS.GSpi, false);
      CalcPointJacobian6D(model, Q, CS.body_s[c], CS.X_s[c].r, CS.GSsi, false);
      CS.GSJ = CS.GSsi - CS.GSpi;

      // Compute position and rotation matrix from predecessor body to base.
      pos_p = CalcBodyToBaseCoordinates (model, Q, CS.body_p[c], CS.X_p[c].r
          , false);
      rot_p = CalcBodyWorldOrientation (model, Q, CS.body_p[c]
          , false).transpose()* CS.X_p[c].E;
      X_0p = SpatialTransform (rot_p, pos_p);

      // Update variables for optimization check.
      prev_constraint_type = ConstraintSet::LoopConstraint;
      prev_body_id_1 = CS.body_p[c];
      prev_body_id_2 = CS.body_s[c];
      prev_body_X_1 = CS.X_p[c];
      prev_body_X_2 = CS.X_s[c];
    }

    // Express the constraint axis in the base frame.
    axis = X_0p.apply(CS.constraintAxis[c]);

    // Compute the constraint Jacobian row.
    G.block(c, 0, 1, model.dof_count) = axis.transpose() * CS.GSJ;
  */
  }
}

RBDL_DLLAPI void CalcConstrainedSystemVariables (
    Model &model,
    EDModel &ed_model,
    const VectorNd  &q,
    const MatrixNd  &q_dirs,
    const VectorNd  &qdot,
    const MatrixNd  &qdot_dirs,
    const VectorNd  &tau,
    const MatrixNd  &tau_dirs,
    ConstraintSet   &CS,
    EDConstraintSet &ed_CS
) {
  const unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());
  ed_model.resize_directions(ndirs);
  ed_CS.resize_directions(ndirs);

  // Compute C
  RigidBodyDynamics::ED::NonlinearEffects (
    model, ed_model, q, q_dirs, qdot, qdot_dirs, CS.C, ed_CS.C
  );

  // Compute H
  assert (CS.H.cols() == model.dof_count && CS.H.rows() == model.dof_count);
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ed_CS.H[idir].setZero();
  }
  CS.H.setZero();
  RigidBodyDynamics::ED::CompositeRigidBodyAlgorithm (
    model, ed_model, q, q_dirs, CS.H, ed_CS.H, false
  );

  // Compute G
  // We have to update model.X_base as they are not automatically computed
  // by NonlinearEffects()
  UpdateKinematicsCustom(
    model, ed_model,
    &q, &q_dirs,
    NULL, NULL,
    NULL, NULL
  );

  // TODO: !!! This call to UpdateKinematicsCustom writes an invalid value
  //           into model.v_J to an invalid value because inside, jcalc is
  //           without qdot and qdot_dirs !!! CHANGE!
  // TODO replace with exact code, saves computation of model.X_lambda
  // NOTE we have to compute ed_model.X_lambda on top
  // for (unsigned int i = 1; i < model.mBodies.size(); i++) {
  //   unsigned int lambda = model.lambda[i];
  //   for (unsigned idir = 0; idir < ndirs; idir++) {
  //     // derivative evaluation
  //     ed_model.X_base[i][idir] = SpatialTransform (
  //       // E * XT.E,
  //       ed_model.X_lambda[i][idir].E * model.X_base[lambda].E
  //       + model.X_lambda[i].E * ed_model.X_base[lambda][idir].E,
  //       // XT.r + XT.E.transpose() * r
  //       ed_model.X_base[lambda][idir].r
  //       + ed_model.X_base[lambda][idir].E.transpose() * model.X_lambda[i].r
  //       + model.X_base[lambda].E.transpose() * ed_model.X_lambda[i][idir].r
  //     );
  //   }
  //   // nominal evaluation
  //   model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
  // }

  RigidBodyDynamics::ED::CalcConstraintsJacobian(
    model, ed_model, q, q_dirs, CS, ed_CS, CS.G, ed_CS.G, false
  );

  /*
  // Compute position error for Baumgarte Stabilization.
  CalcConstraintsPositionError (model, ed_model, q, q_dirs,
                                CS, ed_CS, CS.err, ed_CS.err, false);

  // Compute velocity error for Baugarte stabilization.
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ed_CS.errd.col(idir) = ed_CS.G[idir] * qdot + CS.G * qdot_dirs.col(idir);
  }
  CS.errd = CS.G * qdot;
  */

  // Compute gamma
  unsigned int prev_body_id = 0;
  MatrixNd ed_prev_body_point = MatrixNd::Zero(3, ndirs);
  MatrixNd ed_gamma_i = MatrixNd::Zero(3, ndirs);
  Vector3d prev_body_point = Vector3d::Zero();
  Vector3d gamma_i = Vector3d::Zero();

  ed_CS.QDDot_0.setZero();
  CS.QDDot_0.setZero();
  UpdateKinematicsCustom(
      model, ed_model,   // changed from
      0, &q_dirs, // &q, &q_dirs,       // NULL, &q_dirs,
      0, 0, // &qdot, &qdot_dirs, // NULL, NULL    to get valid values for v_J in model
      &CS.QDDot_0, &ed_CS.QDDot_0 // => this should be improved as it is slow
    );

  CS.gamma = CS.acceleration;
  for (unsigned int i = 0; i < CS.contactConstraintIndices.size(); i++) {
    unsigned const c = CS.contactConstraintIndices[i];

    // only compute point accelerations when necessary
    if (prev_body_id != CS.body[c] || prev_body_point != CS.point[c]) {
      gamma_i = RigidBodyDynamics::ED::CalcPointAcceleration(
            model, ed_model, q, q_dirs, qdot, qdot_dirs,
            CS.QDDot_0, ed_CS.QDDot_0, CS.body[c], CS.point[c], ed_gamma_i,
            false);
      // gamma_i = CalcPointAcceleration (
      //   model, Q, QDot, CS.QDDot_0, CS.body[c], CS.point[c], false
      // );
      prev_body_id = CS.body[c];
      prev_body_point = CS.point[c];
    }

    // we also substract ContactData[c].acceleration such that the contact
    // point will have the desired acceleration
     ed_CS.gamma.row(c).leftCols(ndirs) = -CS.normal[c].transpose() * ed_gamma_i;
     CS.gamma[c] -= CS.normal[c].dot(gamma_i);
  }



  // LOOP CONSTRAINTS NOT SUPPORTED
  /*
  for (unsigned int i = 0; i < CS.loopConstraintIndices.size(); i++) {
    const unsigned int c = CS.loopConstraintIndices[i];

    // Variables used for computations.
    MatrixNd ed_pos_p(3, ndirs);
    vector<Matrix3d> ed_rot_p(ndirs);
    vector<SpatialVector> ed_vel_p(ndirs);
    vector<SpatialVector> ed_vel_s(ndirs);
    vector<SpatialVector> ed_axis(ndirs);

    Vector3d pos_p;
    Matrix3d rot_p;
    SpatialVector vel_p;
    SpatialVector vel_s;
    SpatialVector axis;
    unsigned int id_p;
    unsigned int id_s;

    // Force recomputation.
    prev_body_id = 0;

    // Express the constraint axis in the base frame.
    // nominal + derivative
    pos_p = CalcBodyToBaseCoordinates(model, ed_model, q, q_dirs, CS.body_p[c],
                                      CS.X_p[c].r, ed_pos_p, false);
    // nominal + derivative
    rot_p = CalcBodyWorldOrientation(model, ed_model, q, q_dirs, CS.body_p[c],
                                     ed_rot_p, false);
    // derivative
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ed_rot_p[idir] = ed_rot_p[idir].transpose() * CS.X_p[c].E;
    }
    // nominal
    rot_p = rot_p.transpose() * CS.X_p[c].E;

    // derivative + nominal
    SpatialTransform st(rot_p, pos_p);
    vector<SpatialTransform> st_dirs(ndirs);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      st_dirs[idir] = SpatialTransform(ed_rot_p[idir], ed_pos_p.col(idir));
    }
    applySTSV(ndirs,
              st, st_dirs,
              CS.constraintAxis[c],
              axis, ed_axis);

    // Compute the spatial velocities of the two constrained bodies.
    vel_p = CalcPointVelocity6D (model, ed_model, q, q_dirs, qdot, qdot_dirs,
                                 CS.body_p[c], CS.X_p[c].r,
                                 ed_vel_p, false);

    vel_s = CalcPointVelocity6D (model, ed_model, q, q_dirs, qdot, qdot_dirs,
                                 CS.body_s[c], CS.X_s[c].r,
                                 ed_vel_s, false);

    // Check if the bodies involved in the constraint are fixed. If yes, find
    // their movable parent to access the right value in the a vector.
    // This is needed because we access the model.a vector directly later.
    id_p = GetMovableBodyId (model, CS.body_p[c]);
    id_s = GetMovableBodyId (model, CS.body_s[c]);

    // Problem here if one of the bodies is fixed...
    // Compute the value of gamma.

    SpatialVector temp_accel =
        model.a[id_s] - model.a[id_p] + Math::crossm(vel_s, vel_p);

    // derivative
    for (unsigned idir = 0; idir < ndirs; idir++) {

      SpatialVector ed_temp_accel =
          ed_model.a[id_s][idir]
          - ed_model.a[id_p][idir]
          + AD::crossm(vel_s, ed_vel_s[idir],
                       vel_p, ed_vel_p[idir]);

      double ed_axis_temp_accel = ed_axis[idir].transpose() * temp_accel;

      ed_CS.gamma(c, idir)
        = - axis.transpose() * ed_temp_accel
          - ed_axis_temp_accel
          - 2. * CS.T_stab_inv[c] * ed_CS.errd(c, idir)
          - CS.T_stab_inv[c] * CS.T_stab_inv[c] * ed_CS.err(c, idir);
    }

    // nominal
    CS.gamma[c]
      // Right hand side term.
      = - axis.transpose() * temp_accel
      // Baumgarte stabilization term.
      - 2. * CS.T_stab_inv[c] * CS.errd[c]
      - CS.T_stab_inv[c] * CS.T_stab_inv[c] * CS.err[c];
  }
  */
}


RBDL_DLLAPI void ForwardDynamicsConstraintsDirect (
    Model   &model,
    EDModel &ed_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs,
    const VectorNd &tau,
    const MatrixNd &tau_dirs,
    ConstraintSet &CS,
    EDConstraintSet &ed_CS,
    VectorNd &qddot,
    MatrixNd &ed_qddot) {
  LOG << "-------- " << __func__ << " --------" << std::endl;
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == static_cast<unsigned>(qdot_dirs.cols()));
  assert(ndirs == static_cast<unsigned>(tau_dirs.cols()));

  ed_model.resize_directions(ndirs);
  ed_CS.resize_directions(ndirs);

  // derivative and nominal code
  CalcConstrainedSystemVariables (
      model, ed_model,
      q, q_dirs,
      qdot, qdot_dirs,
      tau, tau_dirs,
      CS, ed_CS
  );

  // derivative and nominal code
  VectorNd c    = tau - CS.C;
  MatrixNd ed_c(c.rows(), ndirs);
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ed_c.col(idir) = tau_dirs.col(idir) - ed_CS.C.col(idir);
  }
  RigidBodyDynamics::ED::SolveConstrainedSystemDirect (
    CS.H, ed_CS.H, CS.G, ed_CS.G, c, ed_c,
    CS.gamma, ed_CS.gamma, CS.A, ed_CS.A, CS.b, ed_CS.b,
    CS.x, ed_CS.x, CS.linear_solver, ndirs
  );

  // Copy back QDDot
  // derivative code
  ed_qddot = ed_CS.x.block(0, 0, model.dof_count, ndirs);
  // nominal code
  qddot = CS.x.segment(0, model.dof_count);

  // Copy back contact forces
  // derivative code
  ed_CS.force = -ed_CS.x.block(model.dof_count, 0, CS.size(), ndirs);
  // nominal code
  CS.force    = -CS.x.segment(model.dof_count, CS.size());
}

// -----------------------------------------------------------------------------
} // namespace ED
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
