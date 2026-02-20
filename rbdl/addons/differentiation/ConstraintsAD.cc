#include "rbdl/Constraints.h"
#include "ConstraintsAD.h"
#include "DynamicsAD.h"
#include "rbdl/SpatialAlgebraOperators.h"
#include "SpatialAlgebraOperatorsAD.h"
#include "rbdl_mathutilsAD.h"

#include <iomanip>

using namespace std;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;

ADConstraintSet::ADConstraintSet(const ConstraintSet & CS_, int _dof_count)
  : ndirs (0)
  , dof_count (0)
  , CS (0)
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

void ADConstraintSet::resize_directions(const unsigned int& requested_ndirs) {
  // std::cout << "in " << __func__ << std::endl;
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
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void ComputeConstraintImpulsesDirect (
    Model &model,
    ADModel &ad_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot_minus,
    const MatrixNd &qdot_minus_dirs,
    ConstraintSet &CS,
    ADConstraintSet &ad_CS,
    VectorNd &qdot_plus,
    MatrixNd &ad_qdot_plus) {
  int ndirs = q_dirs.cols();
  assert(ndirs == qdot_minus_dirs.cols());
  assert(ndirs == ad_qdot_plus.cols());

  // Compute H
  UpdateKinematicsCustom(model, ad_model, &q, &q_dirs, 0, 0, 0, 0);
  CompositeRigidBodyAlgorithm(model, ad_model, q, q_dirs, CS.H, ad_CS.H, false);

  // Compute G
  CalcConstraintsJacobian (model, ad_model, q, q_dirs,
                       CS, ad_CS, CS.G, ad_CS.G, false);

  VectorNd c = CS.H * qdot_minus;
  MatrixNd c_ad(CS.H.rows(), ndirs);
  for (int i = 0; i < ndirs; i++) {
    c_ad.block(0, i, CS.H.rows(), 1) = CS.H * qdot_minus_dirs.col(i) + ad_CS.H[i] * qdot_minus;
  }

  SolveConstrainedSystemDirect (CS.H, ad_CS.H, CS.G, ad_CS.G, c, c_ad, CS.v_plus,
                            ad_CS.v_plus, CS.A, ad_CS.A, CS.b, ad_CS.b,
                            CS.x, ad_CS.x, CS.linear_solver, ndirs);

  // Copy back QDotPlus
  // derivative evaluation
  ad_qdot_plus = ad_CS.x.block(0, 0, model.dof_count, ndirs);
  // nominal evaluation
  qdot_plus = CS.x.segment(0, model.dof_count);

  // Copy back constraint impulses
  ad_CS.impulse = ad_CS.x.block(model.dof_count, 0, CS.size(), ndirs);
  CS.impulse = CS.x.segment(model.dof_count, CS.size());
}

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
  LOG << "x = " << std::endl << x << std::endl;
}

RBDL_DLLAPI void CalcConstraintsPositionError (
    Model &model,
    ADModel &ad_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    ConstraintSet &cs,
    ADConstraintSet &ad_cs,
    VectorNd &err,
    MatrixNd &ad_err,
    bool update_kinematics) {
  unsigned const ndirs = q_dirs.cols();
  assert(static_cast<unsigned>(err.size()) == cs.size());
  assert(ndirs <= ad_err.cols());

  if(update_kinematics) {
    UpdateKinematicsCustom(model, ad_model,
                           &q, &q_dirs,
                           NULL, NULL,
                           NULL, NULL);
  }

  for (unsigned int i = 0; i < cs.contactConstraintIndices.size(); i++) {
    unsigned const c = cs.contactConstraintIndices[i];
    ad_err.row(c).setZero();
    err[c] = 0.;
  }

  for (unsigned int i = 0; i < cs.loopConstraintIndices.size(); i++) {
    unsigned const lci = cs.loopConstraintIndices[i];

    MatrixNd ad_pos_p(3, ndirs);
    MatrixNd ad_pos_s(3, ndirs);
    vector<Matrix3d> ad_rot_p(ndirs);
    vector<Matrix3d> ad_rot_s(ndirs);
    vector<Matrix3d> ad_rot_ps(ndirs);
    vector<SpatialVector> ad_d(ndirs);

    // Variables used for computations.
    Vector3d pos_p;
    Vector3d pos_s;
    Matrix3d rot_p;
    Matrix3d rot_s;
    Matrix3d rot_ps;
    SpatialVector d;

    // Constraints computed in the predecessor body frame.

    // Compute the orientation of the two constraint frames.
    rot_p = CalcBodyWorldOrientation (model, ad_model, q, q_dirs,
                                      cs.body_p[lci], ad_rot_p, false);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_rot_p[idir] = ad_rot_p[idir].transpose() * cs.X_p[lci].E;
    }
    rot_p = rot_p.transpose() * cs.X_p[lci].E;

    rot_s = CalcBodyWorldOrientation (model, ad_model, q, q_dirs,
                                      cs.body_s[lci], ad_rot_s, false);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_rot_s[idir] = ad_rot_s[idir].transpose() * cs.X_p[lci].E;
    }
    rot_s = rot_s.transpose() * cs.X_s[lci].E;

    // Compute the orientation from the predecessor to the successor frame.
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_rot_ps[idir] =
          ad_rot_p[idir].transpose() * rot_s
          + rot_p.transpose() * ad_rot_s[idir];
    }
    rot_ps = rot_p.transpose() * rot_s;

    // Compute the position of the two contact points.
    pos_p = CalcBodyToBaseCoordinates(model, ad_model, q, q_dirs,
                                      cs.body_p[lci], cs.X_p[lci].r,
                                      ad_pos_p, false);
    pos_s = CalcBodyToBaseCoordinates(model, ad_model, q, q_dirs,
                                      cs.body_s[lci], cs.X_s[lci].r,
                                      ad_pos_s, false);

    // The first three elemenets represent the rotation error.
    // This formulation is equivalent to u * sin(theta), where u and theta are
    // the angle-axis of rotation from the predecessor to the successor frame.
    // These quantities are expressed in the predecessor frame.
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_d[idir](0) = -0.5 * (ad_rot_ps[idir](1,2) - ad_rot_ps[idir](2,1));
      ad_d[idir](1) = -0.5 * (ad_rot_ps[idir](2,0) - ad_rot_ps[idir](0,2));
      ad_d[idir](2) = -0.5 * (ad_rot_ps[idir](0,1) - ad_rot_ps[idir](1,0));
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
      ad_d[idir].segment<3>(3) =
          ad_rot_p[idir].transpose() * (pos_s - pos_p)
          + rot_p.transpose() * (ad_pos_s.col(idir) - ad_pos_p.col(idir));
    }
    d.segment<3>(3) = rot_p.transpose() * (pos_s - pos_p);

    // Project the error on the constraint axis to find the actual error.
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_err.col(idir) = cs.constraintAxis[lci].transpose() * ad_d[idir];
    }
    err[lci] = cs.constraintAxis[lci].transpose() * d;
  }
}

RBDL_DLLAPI void CalcConstraintsJacobian(
    Model &model,
    ADModel &ad_model,
    const VectorNd   &q,
    const MatrixNd   &q_dirs,
    ConstraintSet    &CS,
    ADConstraintSet  &ad_CS,
    MatrixNd         &G,
    vector<MatrixNd> &G_dirs,
    bool update_kinematics
) {
  const unsigned int ndirs = q_dirs.cols();
  assert(ndirs <= G_dirs.size());

  // resize if required
  ad_model.resize_directions(ndirs);
  ad_CS.resize_directions(ndirs);

  if (update_kinematics) {
    // nominal + derivative evaluation
    UpdateKinematicsCustom (model, ad_model, &q, &q_dirs, 0, 0, 0, 0);
  }

  // variables to check whether we need to recompute G.
  ConstraintSet::ConstraintType prev_constraint_type
    = ConstraintSet::ConstraintTypeLast;
  unsigned prev_body_id_1 = 0;
  unsigned prev_body_id_2 = 0;

  vector<SpatialTransform> prev_body_X_1_dirs(ndirs);
  vector<SpatialTransform> prev_body_X_2_dirs(ndirs);
  SpatialTransform prev_body_X_1;
  SpatialTransform prev_body_X_2;

  for (unsigned int i = 0; i < CS.contactConstraintIndices.size(); i++) {
    const unsigned int c = CS.contactConstraintIndices[i];

    // only compute the matrix Gi if actually needed
    if (prev_constraint_type != CS.constraintType[c]
        || prev_body_id_1 != CS.body[c]
        || prev_body_X_1.r != CS.point[c]) {

      // Compute the jacobian for the point.
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_CS.Gi[idir].setZero();
      }
      CS.Gi.setZero();
      CalcPointJacobian (model, ad_model, q, q_dirs,
                         CS.body[c], CS.point[c],
                         CS.Gi, ad_CS.Gi,
                         false);

      prev_constraint_type = ConstraintSet::ContactConstraint;

      // Update variables for optimization check.
      prev_body_id_1 = CS.body[c];
      for (unsigned idir = 0; idir < ndirs; idir++) {
        prev_body_X_1_dirs[idir].setZero();
      }
      prev_body_X_1 = RigidBodyDynamics::Xtrans(CS.point[c]);
    }

    for(unsigned int j = 0; j < model.dof_count; j++) {
      MatrixNd gaxis_dirs(3, ndirs);
      for (unsigned idir = 0; idir < ndirs; idir++) {
        gaxis_dirs.col(idir) = ad_CS.Gi[idir].col(j);
      }
      Vector3d gaxis (CS.Gi(0,j), CS.Gi(1,j), CS.Gi(2,j));
      for (unsigned idir = 0; idir < ndirs; idir++) {
        G_dirs[idir](c,j) = gaxis_dirs.col(idir).transpose() * CS.normal[c];
      }
      G(c,j) = gaxis.transpose() * CS.normal[c];
    }
  }

  // Variables used for computations.
  MatrixNd ad_normal(3, ndirs);
  vector<SpatialVector> ad_axis(ndirs);
  MatrixNd ad_pos_p(3, ndirs);
  vector<Matrix3d> ad_rot_p(ndirs);
  vector<SpatialTransform> ad_X_0p(ndirs);

  Vector3d normal;
  SpatialVector axis;
  Vector3d pos_p;
  Matrix3d rot_p;
  SpatialTransform X_0p;

  for (unsigned i = 0; i < CS.loopConstraintIndices.size(); i++) {
    unsigned const c = CS.loopConstraintIndices[i];

    // Only recompute variables if necessary.
    if( prev_body_id_1 != CS.body_p[c]
        || prev_body_id_2 != CS.body_s[c]
        || prev_body_X_1.r != CS.X_p[c].r
        || prev_body_X_2.r != CS.X_s[c].r
        || prev_body_X_1.E != CS.X_p[c].E
        || prev_body_X_2.E != CS.X_s[c].E) {

      // Compute the 6D jacobians of the two contact points.
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_CS.GSpi[idir].setZero();
      }
      CS.GSpi.setZero();

      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_CS.GSsi[idir].setZero();
      }
      CS.GSsi.setZero();

      CalcPointJacobian6D(model, ad_model, q, q_dirs, CS.body_p[c],
                          CS.X_p[c].r, CS.GSpi, ad_CS.GSpi, false);

      CalcPointJacobian6D(model, ad_model, q, q_dirs, CS.body_s[c],
                          CS.X_s[c].r, CS.GSsi, ad_CS.GSsi, false);

      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_CS.GSJ[idir] = ad_CS.GSsi[idir] - ad_CS.GSpi[idir];
      }
      CS.GSJ = CS.GSsi - CS.GSpi;

      // Compute position and rotation matrix from predecessor body to base.

      pos_p = CalcBodyToBaseCoordinates(model, ad_model, q, q_dirs,
                                        CS.body_p[c], CS.X_p[c].r,
                                        ad_pos_p, false);

      rot_p = CalcBodyWorldOrientation(model, ad_model, q, q_dirs,
                                       CS.body_p[c], ad_rot_p, false);
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_rot_p[idir] = ad_rot_p[idir].transpose() * CS.X_p[c].E;
      }
      rot_p = rot_p.transpose() * CS.X_p[c].E;


      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_X_0p[idir] = SpatialTransform(ad_rot_p[idir], ad_pos_p.col(idir));
      }
      X_0p = SpatialTransform (rot_p, pos_p);

      // Update variables for optimization check.
      prev_constraint_type = ConstraintSet::LoopConstraint;
      prev_body_id_1 = CS.body_p[c];
      prev_body_id_2 = CS.body_s[c];
      prev_body_X_1 = CS.X_p[c];
      prev_body_X_2 = CS.X_s[c];
    }

    // Express the constraint axis in the base frame.
    applySTSV(ndirs,
              X_0p, ad_X_0p,
              CS.constraintAxis[c],
              axis, ad_axis);

    // Compute the constraint Jacobian row.
    for (unsigned idir = 0; idir < ndirs; idir++) {
      G_dirs[idir].row(c) =
          axis.transpose() * ad_CS.GSJ[idir]
          + ad_axis[idir].transpose() * CS.GSJ;
    }
    G.row(c) = axis.transpose() * CS.GSJ;
  }
}

RBDL_DLLAPI void CalcConstrainedSystemVariables (
    Model &model,
    ADModel &ad_model,
    const VectorNd  &q,
    const MatrixNd  &q_dirs,
    const VectorNd  &qdot,
    const MatrixNd  &qdot_dirs,
    const VectorNd  &tau,
    const MatrixNd  &tau_dirs,
    ConstraintSet   &CS,
    ADConstraintSet &ad_CS) {
  unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == tau_dirs.cols());
  ad_model.resize_directions(ndirs);
  ad_CS.resize_directions(ndirs);

  // Compute C
  NonlinearEffects (model, ad_model, q, q_dirs, qdot, qdot_dirs, CS.C, ad_CS.C);
  assert (CS.H.cols() == model.dof_count && CS.H.rows() == model.dof_count);

  // Compute H
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_CS.H[idir].setZero();
  }
  CS.H.setZero();
  CompositeRigidBodyAlgorithm (model, ad_model, q, q_dirs, CS.H, ad_CS.H, false);

  // Compute G
  // We have to update model.X_base as they are not automatically computed
  // by NonlinearEffects()
  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    unsigned int lambda = model.lambda[i];
    mulSTST (ndirs,
             model.X_lambda[i], ad_model.X_lambda[i],
             model.X_base[lambda], ad_model.X_base[lambda],
             model.X_base[i], ad_model.X_base[i]);
  }
  CalcConstraintsJacobian(model, ad_model, q, q_dirs, CS, ad_CS, CS.G, ad_CS.G, false);

  // Compute position error for Baumgarte Stabilization.
  CalcConstraintsPositionError (model, ad_model, q, q_dirs,
                                CS, ad_CS, CS.err, ad_CS.err, false);

  // Compute velocity error for Baugarte stabilization.
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_CS.errd.col(idir) = ad_CS.G[idir] * qdot + CS.G * qdot_dirs.col(idir);
  }
  CS.errd = CS.G * qdot;

  // Compute gamma
  unsigned int prev_body_id = 0;
  MatrixNd ad_prev_body_point = MatrixNd::Zero(3, ndirs);
  MatrixNd ad_gamma_i = MatrixNd::Zero(3, ndirs);
  Vector3d prev_body_point = Vector3d::Zero();
  Vector3d gamma_i = Vector3d::Zero();

  ad_CS.QDDot_0.setZero();
  CS.QDDot_0.setZero();
  UpdateKinematicsCustom(model, ad_model,
                         NULL, NULL,
                         NULL, NULL,
                         &CS.QDDot_0, &ad_CS.QDDot_0);

  for (unsigned int i = 0; i < CS.contactConstraintIndices.size(); i++) {
    unsigned const c = CS.contactConstraintIndices[i];

    // only compute point accelerations when necessary
    if (prev_body_id != CS.body[c] || prev_body_point != CS.point[c]) {
      gamma_i = CalcPointAcceleration(
            model, ad_model, q, q_dirs, qdot, qdot_dirs,
            CS.QDDot_0, ad_CS.QDDot_0, CS.body[c], CS.point[c], ad_gamma_i,
            false);
      prev_body_id = CS.body[c];
      prev_body_point = CS.point[c];
    }

    // we also substract ContactData[c].acceleration such that the contact
    // point will have the desired acceleration
    ad_CS.gamma.row(c).segment(0, ndirs) =
        -CS.normal[c].transpose() * ad_gamma_i;
    CS.gamma[c] = CS.acceleration[c] - CS.normal[c].dot(gamma_i);
  }

  for (unsigned int i = 0; i < CS.loopConstraintIndices.size(); i++) {
    const unsigned int c = CS.loopConstraintIndices[i];

    // Variables used for computations.
    MatrixNd ad_pos_p(3, ndirs);
    vector<Matrix3d> ad_rot_p(ndirs);
    vector<SpatialVector> ad_vel_p(ndirs);
    vector<SpatialVector> ad_vel_s(ndirs);
    vector<SpatialVector> ad_axis(ndirs);

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
    pos_p = CalcBodyToBaseCoordinates(model, ad_model, q, q_dirs, CS.body_p[c],
                                      CS.X_p[c].r, ad_pos_p, false);
    // nominal + derivative
    rot_p = CalcBodyWorldOrientation(model, ad_model, q, q_dirs, CS.body_p[c],
                                     ad_rot_p, false);
    // derivative
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_rot_p[idir] = ad_rot_p[idir].transpose() * CS.X_p[c].E;
    }
    // nominal
    rot_p = rot_p.transpose() * CS.X_p[c].E;

    // derivative + nominal
    SpatialTransform st(rot_p, pos_p);
    vector<SpatialTransform> st_dirs(ndirs);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      st_dirs[idir] = SpatialTransform(ad_rot_p[idir], ad_pos_p.col(idir));
    }
    applySTSV(ndirs,
              st, st_dirs,
              CS.constraintAxis[c],
              axis, ad_axis);

    // Compute the spatial velocities of the two constrained bodies.
    vel_p = CalcPointVelocity6D (model, ad_model, q, q_dirs, qdot, qdot_dirs,
                                 CS.body_p[c], CS.X_p[c].r,
                                 ad_vel_p, false);

    vel_s = CalcPointVelocity6D (model, ad_model, q, q_dirs, qdot, qdot_dirs,
                                 CS.body_s[c], CS.X_s[c].r,
                                 ad_vel_s, false);

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

      SpatialVector ad_temp_accel =
          ad_model.a[id_s][idir]
          - ad_model.a[id_p][idir]
          + AD::crossm(vel_s, ad_vel_s[idir],
                       vel_p, ad_vel_p[idir]);

      double ad_axis_temp_accel = ad_axis[idir].transpose() * temp_accel;

      ad_CS.gamma(c, idir)
        = - axis.transpose() * ad_temp_accel
          - ad_axis_temp_accel
          - 2. * CS.T_stab_inv[c] * ad_CS.errd(c, idir)
          - CS.T_stab_inv[c] * CS.T_stab_inv[c] * ad_CS.err(c, idir);
    }

    // nominal
    CS.gamma[c]
      // Right hand side term.
      = - axis.transpose() * temp_accel
      // Baumgarte stabilization term.
      - 2. * CS.T_stab_inv[c] * CS.errd[c]
      - CS.T_stab_inv[c] * CS.T_stab_inv[c] * CS.err[c];
  }
}

RBDL_DLLAPI void ForwardDynamicsConstraintsDirect (
    Model   &model,
    ADModel &ad_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs,
    const VectorNd &tau,
    const MatrixNd &tau_dirs,
    ConstraintSet &CS,
    ADConstraintSet &ad_CS,
    VectorNd &qddot,
    MatrixNd &ad_qddot) {
  LOG << "-------- " << __func__ << " --------" << std::endl;
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == static_cast<unsigned>(qdot_dirs.cols()));
  assert(ndirs == static_cast<unsigned>(tau_dirs.cols()));

  ad_model.resize_directions(ndirs);
  ad_CS.resize_directions(ndirs);

  // derivative and nominal code
  CalcConstrainedSystemVariables (
      model, ad_model,
      q, q_dirs,
      qdot, qdot_dirs,
      tau, tau_dirs,
      CS, ad_CS
  );

  // derivative and nominal code
  VectorNd c    = tau - CS.C;
  MatrixNd ad_c(c.rows(), ndirs);
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_c.col(idir) = tau_dirs.col(idir) - ad_CS.C.col(idir);
  }
  SolveConstrainedSystemDirect (CS.H, ad_CS.H, CS.G, ad_CS.G, c, ad_c,
                            CS.gamma, ad_CS.gamma, CS.A, ad_CS.A, CS.b, ad_CS.b,
                            CS.x, ad_CS.x, CS.linear_solver, ndirs);

  // Copy back QDDot
  // derivative code
  ad_qddot = ad_CS.x.block(0, 0, model.dof_count, ndirs);
  // nominal code
  qddot = CS.x.segment(0, model.dof_count);

  // Copy back contact forces
  // derivative code
  ad_CS.force = -ad_CS.x.block(model.dof_count, 0, CS.size(), ndirs);
  // nominal code
  CS.force    = -CS.x.segment(model.dof_count, CS.size());
}

RBDL_DLLAPI
void ForwardDynamicsAccelerationDeltas (
    Model &model,
    ADModel &ad_model,
    ConstraintSet &CS,
    ADConstraintSet &ad_CS,
    VectorNd &QDDot_t,
    MatrixNd &ad_QDDot_t,
    const unsigned int body_id,
    const std::vector<SpatialVector> &f_t,
    const std::vector<MatrixNd> &f_t_dirs
) {
  LOG << "-------- " << __func__ << " ------" << std::endl;

  const unsigned ndirs = ad_QDDot_t.cols();
  ad_model.resize_directions(ndirs);
  ad_CS.resize_directions(ndirs);

  // TODO reset all values (debug)
  for (unsigned int i = 0; i < model.mBodies.size(); i++) {
    ad_CS.d_pA[i].setZero();
    ad_CS.d_a[i].setZero();
    ad_CS.d_u.setZero();
    ad_CS.d_multdof3_u[i].setZero();

    CS.d_pA[i].setZero();
    CS.d_a[i].setZero();
    CS.d_u[i] = 0.;
    CS.d_multdof3_u[i].setZero();
  }
  for(unsigned int i=0; i<model.mCustomJoints.size();i++){
    // TODO expose in ADModel
    model.mCustomJoints[i]->d_u.setZero();
  }

  for (unsigned int i = body_id; i > 0; i--) {
    if (i == body_id) {
      // derivative evaluation
      applyAdjointSTSV (
        ndirs,
        model.X_base[i], ad_model.X_base[i],
        -f_t[i], -f_t_dirs[i],
        CS.d_pA[i], ad_CS.d_pA[i]
      );
      // nominal evaluation
      // CS.d_pA[i] = -model.X_base[i].applyAdjoint(f_t[i]);
    }

    if (model.mJoints[i].mDoFCount == 3
        && model.mJoints[i].mJointType != JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__
           << ": multi-dof joints not supported." << endl;
      abort();

      // CS.d_multdof3_u[i] = - model.multdof3_S[i].transpose() * (CS.d_pA[i]);

      // unsigned int lambda = model.lambda[i];
      // if (lambda != 0) {
      //   CS.d_pA[lambda] =   CS.d_pA[lambda]
      //     + model.X_lambda[i].applyTranspose (
      //         CS.d_pA[i] + (model.multdof3_U[i]
      //           * model.multdof3_Dinv[i]
      //           * CS.d_multdof3_u[i]));
      // }
    } else if(model.mJoints[i].mDoFCount == 1
        && model.mJoints[i].mJointType != JointTypeCustom) {

      // derivative evaluation
      for (unsigned int idir = 0; idir < ndirs; ++idir)
      {
        ad_CS.d_u(i, idir) = - model.S[i].dot(ad_CS.d_pA[i].col(idir));
      }
      // nominal evaluation
      CS.d_u[i] = - model.S[i].dot(CS.d_pA[i]);
      unsigned int lambda = model.lambda[i];

      if (lambda != 0) {
        // derivative evaluation
        MatrixNd temp_dirs = MatrixNd::Zero(6, ndirs);
        for (unsigned int idir = 0; idir < ndirs; ++idir)
        {
          temp_dirs.col(idir)
            = ad_CS.d_pA[i].col(idir)
            + ad_model.U[i][idir] * CS.d_u[i] / model.d[i]
            + model.U[i] * ad_CS.d_u(i, idir) / model.d[i]
            - model.U[i] * CS.d_u[i] * ad_model.d(i, idir) / (model.d[i] * model.d[i]);
          ;
        }
        addApplyTransposeSTSV(
          ndirs,
          1.0,
          // SpatialTransform X, X_dirs
          model.X_lambda[i], ad_model.X_lambda[i],
          // SpatialVector v, v_dirs
          CS.d_pA[i] + model.U[i] * CS.d_u[i] / model.d[i], temp_dirs,
          // res, res_dirs
          CS.d_pA[lambda], ad_CS.d_pA[lambda]
        );
        // nominal evaluation
        // CS.d_pA[lambda] += model.X_lambda[i].applyTranspose (
        //       CS.d_pA[i] + model.U[i] * CS.d_u[i] / model.d[i]
        // );
      }
    } else if (model.mJoints[i].mJointType == JointTypeCustom){
      cerr << __FILE__ << " " << __LINE__
           << ": Custom joints not supported." << endl;
      abort();

      // unsigned int kI     = model.mJoints[i].custom_joint_index;
      // // unsigned int dofI   = model.mCustomJoints[kI]->mDoFCount;
      // //CS.
      // model.mCustomJoints[kI]->d_u =
      //   - model.mCustomJoints[kI]->S.transpose() * (CS.d_pA[i]);
      // unsigned int lambda = model.lambda[i];
      // if (lambda != 0) {
      //   CS.d_pA[lambda] =
      //     CS.d_pA[lambda]
      //     + model.X_lambda[i].applyTranspose (
      //         CS.d_pA[i] + (   model.mCustomJoints[kI]->U
      //           * model.mCustomJoints[kI]->Dinv
      //           * model.mCustomJoints[kI]->d_u)
      //         );
      // }
    }
  }

  for (unsigned int i = 0; i < f_t.size(); i++) {
    LOG << "f_t[" << i << "] = " << f_t[i].transpose() << std::endl;
  }

  for (unsigned int i = 0; i < model.mBodies.size(); i++) {
    LOG << "i = " << i << ": d_pA[i] " << CS.d_pA[i].transpose() << std::endl;
  }
  for (unsigned int i = 0; i < model.mBodies.size(); i++) {
    LOG << "i = " << i << ": d_u[i] = " << CS.d_u[i] << std::endl;
  }

  // derivative evaluation
  ad_QDDot_t.row(0).setZero();
  for (unsigned int idir = 0; idir < ndirs; ++idir)
  {
    ad_CS.d_a[0].col(idir) = ad_model.a[0][idir];
  }
  std::vector<SpatialVector> ad_Xa
    = std::vector<SpatialVector> (ndirs, SpatialVector::Zero ());
  // nominal evaluation
  QDDot_t[0] = 0.;
  CS.d_a[0] = model.a[0];
  SpatialVector Xa;

  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    unsigned int q_index = model.mJoints[i].q_index;
    unsigned int lambda = model.lambda[i];

    // derivative evaluation
    applySTSV(
      ndirs,
      model.X_lambda[i], ad_model.X_lambda[i],
      CS.d_a[lambda], ad_CS.d_a[lambda],
      Xa, ad_Xa
    );
    // nominal evaluation
    // SpatialVector Xa = model.X_lambda[i].apply(CS.d_a[lambda]);

    if (model.mJoints[i].mDoFCount == 3
        && model.mJoints[i].mJointType != JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__
           << ": multi-dof joints not supported." << endl;
      abort();

      // Vector3d qdd_temp = model.multdof3_Dinv[i]
      //   * (CS.d_multdof3_u[i] - model.multdof3_U[i].transpose() * Xa);

      // QDDot_t[q_index] = qdd_temp[0];
      // QDDot_t[q_index + 1] = qdd_temp[1];
      // QDDot_t[q_index + 2] = qdd_temp[2];
      // model.a[i] = model.a[i] + model.multdof3_S[i] * qdd_temp;
      // CS.d_a[i] = Xa + model.multdof3_S[i] * qdd_temp;
    } else if (model.mJoints[i].mDoFCount == 1
        && model.mJoints[i].mJointType != JointTypeCustom){
      // derivative evaluation
      for (unsigned int idir = 0; idir < ndirs; ++idir)
      {
        ad_QDDot_t(q_index, idir)
          = (
            ad_CS.d_u(i, idir) - ad_model.U[i][idir].dot(Xa) - model.U[i].dot(ad_Xa[idir])
            ) / model.d[i]
          - (CS.d_u[i] - model.U[i].dot(Xa) ) * ad_model.d(i,idir) / (model.d[i] * model.d[i]);
        ad_CS.d_a[i].col(idir) = ad_Xa[idir] + model.S[i] * ad_QDDot_t(q_index, idir);
      }

      // nominal evaluation
      QDDot_t[q_index] = (CS.d_u[i] - model.U[i].dot(Xa) ) / model.d[i];
      CS.d_a[i] = Xa + model.S[i] * QDDot_t[q_index];
    } else if (model.mJoints[i].mJointType == JointTypeCustom){
      cerr << __FILE__ << " " << __LINE__
           << ": Custom joints not supported." << endl;
      abort();

      // unsigned int kI     = model.mJoints[i].custom_joint_index;
      // unsigned int dofI   = model.mCustomJoints[kI]->mDoFCount;
      // VectorNd qdd_temp = VectorNd::Zero(dofI);

      // qdd_temp = model.mCustomJoints[kI]->Dinv
      //   * (model.mCustomJoints[kI]->d_u
      //       - model.mCustomJoints[kI]->U.transpose() * Xa);

      // for(unsigned int z=0; z<dofI;++z){
      //   QDDot_t[q_index+z] = qdd_temp[z];
      // }

      // model.a[i] = model.a[i] + model.mCustomJoints[kI]->S * qdd_temp;
      // CS.d_a[i] = Xa + model.mCustomJoints[kI]->S * qdd_temp;
    }

    LOG << "QDDot_t[" << i - 1 << "] = " << QDDot_t[i - 1] << std::endl;
    LOG << "d_a[i] = " << CS.d_a[i].transpose() << std::endl;
  }
}

RBDL_DLLAPI
void ForwardDynamicsApplyConstraintForces (
    Model &model,
    ADModel &ad_model,
    const VectorNd &Tau,
    const MatrixNd &Tau_dirs,
    ConstraintSet &CS,
    ADConstraintSet &ad_CS,
    VectorNd &QDDot,
    MatrixNd &ad_QDDot
) {
  LOG << "-------- " << __func__ << " --------" << std::endl;
  assert (QDDot.size() == model.dof_count);

  const unsigned int ndirs = Tau_dirs.cols();
  assert(ndirs == static_cast<unsigned>(Tau_dirs.cols()));
  assert(ndirs == static_cast<unsigned>(ad_QDDot.cols()));

  unsigned int i = 0;

  for (i = 1; i < model.mBodies.size(); i++) {
    // derivative code
    for (unsigned int idir = 0; idir < ndirs; idir++) {
      ad_model.IA[i][idir].setZero();
    }
    // nominal code
    model.IA[i] = model.I[i].toMatrix();

    // derivative code
    for (unsigned int idir = 0; idir < ndirs; idir++) {
      ad_model.pA[i][idir] = AD::crossf(
            model.v[i], ad_model.v[i][idir],
            model.I[i] * model.v[i], model.I[i] * ad_model.v[i][idir]);
    }
    // nominal code
    model.pA[i] = Math::crossf(model.v[i], model.I[i] * model.v[i]);

    if (CS.f_ext_constraints[i] != SpatialVector::Zero()) {
      LOG << "External force (" << i << ") = " << model.X_base[i].toMatrixAdjoint() * CS.f_ext_constraints[i] << std::endl;
      addApplyAdjointSTSV(ndirs,
                          -1.0,
                          model.X_base[i], ad_model.X_base[i],
                          CS.f_ext_constraints[i], ad_CS.f_ext_constraints[i],
                          model.pA[i], ad_model.pA[i]);

//      model.pA[i] -= model.X_base[i].toMatrixAdjoint() * CS.f_ext_constraints[i];
    }
  }

  // ClearLogOutput();

  LOG << "--- first loop ---" << std::endl;

  for (i = model.mBodies.size() - 1; i > 0; i--) {
    unsigned int q_index = model.mJoints[i].q_index;

    if (model.mJoints[i].mDoFCount == 3
        && model.mJoints[i].mJointType != JointTypeCustom) {
      std::cerr << "3DoF JOINTS NOT SUPPORTED!" << std::endl;
      abort();
    } else if (model.mJoints[i].mDoFCount == 1
        && model.mJoints[i].mJointType != JointTypeCustom) {
      // derivative code
      for (unsigned int idir = 0; idir < ndirs; idir++) {
        ad_model.u(i, idir) = Tau_dirs(q_index, idir) - model.S[i].dot(ad_model.pA[i][idir]);
      }
      // nominal code
      model.u[i] = Tau[q_index] - model.S[i].dot(model.pA[i]);

      unsigned int lambda = model.lambda[i];
      if (lambda != 0) {
        vector<SpatialMatrix> ad_Ia(ndirs, SpatialMatrix::Zero());
        // derivative evaluation
        for(unsigned int j = 0; j < ndirs; j++) {
          ad_Ia[j] = ad_model.IA[i][j]
            - ad_model.U[i][j] * (model.U[i] / model.d[i]).transpose()
            - model.U[i] * (ad_model.U[i][j] / model.d[i]).transpose()
            - model.U[i] * (model.U[i] * (-ad_model.d(i,j)) / (model.d[i] * model.d[i])).transpose();
        }
        // nominal evaluation
        SpatialMatrix Ia = model.IA[i] - model.U[i] * (model.U[i] / model.d[i]).transpose();

        vector<SpatialVector> ad_pa(ndirs, SpatialVector::Zero());
        // derivative evaluation
        for(unsigned int j = 0; j < ndirs; j++) {
          ad_pa[j] = ad_model.pA[i][j]
            + ad_Ia[j] * model.c[i]
            + Ia * ad_model.c[i][j]
            + ad_model.U[i][j] * model.u[i] / model.d[i]
            + model.U[i] * ad_model.u(i,j) / model.d[i]
            + model.U[i] * model.u[i] * (-ad_model.d(i,j)) / (model.d[i] * model.d[i]);
        }
        // nominal evaluation
        SpatialVector pa = model.pA[i] + Ia * model.c[i] + model.U[i] * model.u[i] / model.d[i];

#ifdef EIGEN_CORE_H
        addSqrFormSTSM_noalias(
              ndirs,
              model.X_lambda[i], ad_model.X_lambda[i],
              Ia, ad_Ia,
              model.IA[lambda], ad_model.IA[lambda]);

        SpatialVector summand;
        vector<SpatialVector> ad_summand(ndirs);
        applyTransposeSTSV(
              ndirs,
              model.X_lambda[i], ad_model.X_lambda[i],
              pa, ad_pa,
              summand, ad_summand);
        for (unsigned idir = 0; idir < ndirs; idir++) {
          ad_model.pA[lambda][idir].noalias() += ad_summand[idir];
        }
        model.pA[lambda].noalias() += summand;
        // original nominal code
        // model.IA[lambda].noalias() += (model.X_lambda[i].toMatrixTranspose() * Ia * model.X_lambda[i].toMatrix());
        // model.pA[lambda].noalias() += model.X_lambda[i].applyTranspose(pa);
#else
        cerr << "Simple math not yet supported." << endl;
        abort();
        // model.IA[lambda] += (model.X_lambda[i].toMatrixTranspose() * Ia * model.X_lambda[i].toMatrix());
        // model.pA[lambda] += model.X_lambda[i].applyTranspose(pa);
#endif
        LOG << "pA[" << lambda << "] = "
          << model.pA[lambda].transpose() << std::endl;
      }
    } else if(model.mJoints[i].mJointType == JointTypeCustom) {
      std::cerr << "Custom joints not supported!" << std::endl;
      abort();
    }
  }

  // derivative code
  for (unsigned int idir = 0; idir < ndirs; idir++) {
    ad_model.a[0][idir].setZero();
  }
  // nominal code
  model.a[0] = SpatialVector (0., 0., 0., -model.gravity[0], -model.gravity[1], -model.gravity[2]);

  for (i = 1; i < model.mBodies.size(); i++) {
    unsigned int q_index = model.mJoints[i].q_index;
    unsigned int lambda = model.lambda[i];

    applySTSV(ndirs,
              model.X_lambda[i], ad_model.X_lambda[i],
              model.a[lambda], ad_model.a[lambda],
              model.a[i], ad_model.a[i]);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_model.a[i][idir] += ad_model.c[i][idir];
    }
    model.a[i] += model.c[i];

//    model.a[i] = X_lambda.apply(model.a[lambda]) + model.c[i];
    LOG << "a'[" << i << "] = " << model.a[i].transpose() << std::endl;

    if (model.mJoints[i].mDoFCount == 3 && model.mJoints[i].mJointType != JointTypeCustom) {
      std::cerr << "3DoF joints not supported!" << std::endl;
      abort();
    } else if (model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom) {
      // derivative evaluation
      for(unsigned idir = 0; idir < ndirs; idir++) {
        ad_QDDot(q_index, idir) =
            -(ad_model.d(i, idir) / (model.d[i] * model.d[i]))
            * (model.u[i] - model.U[i].dot(model.a[i]))
            + (1. / model.d[i])
            * (ad_model.u(i, idir)
               - ad_model.U[i][idir].dot(model.a[i])
               - model.U[i].dot(ad_model.a[i][idir]));
      }
      // nominal evaluation
      QDDot[q_index] = (1./model.d[i]) * (model.u[i] - model.U[i].dot(model.a[i]));

      // derivative evaluation
      for(unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.a[i][idir] = ad_model.a[i][idir] + model.S[i] * ad_QDDot(q_index, idir);
      }
      // nominal evaluation
      model.a[i] = model.a[i] + model.S[i] * QDDot[q_index];
    } else if (model.mJoints[i].mJointType == JointTypeCustom) {
      std::cerr << "Custom joints not supported!" << std::endl;
      abort();
    }
  }

  LOG << "QDDot = " << QDDot.transpose() << std::endl;
}

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
 ) {
  LOG << "-------- " << __func__ << " ------" << std::endl;

  assert (CS.f_ext_constraints.size() == model.mBodies.size());
  assert (CS.QDDot_0.size() == model.dof_count);
  assert (CS.QDDot_t.size() == model.dof_count);
  assert (CS.f_t.size() == CS.size());
  assert (CS.point_accel_0.size() == CS.size());
  assert (static_cast<unsigned>(CS.K.rows()) == CS.size());
  assert (static_cast<unsigned>(CS.K.cols()) == CS.size());
  assert (static_cast<unsigned>(CS.force.size()) == CS.size());
  assert (static_cast<unsigned>(CS.a.size()) == CS.size());

  const unsigned int ndirs = q_dirs.cols();
  assert(ndirs == static_cast<unsigned>(qdot_dirs.cols()));
  assert(ndirs == static_cast<unsigned>(tau_dirs.cols()));

  ad_model.resize_directions(ndirs);
  ad_CS.resize_directions(ndirs);

  Vector3d point_accel_t;

  unsigned int ci = 0;

  // The default acceleration only needs to be computed once
  {
    SUPPRESS_LOGGING;
    // derivative and nominal code
    ForwardDynamics(
      model, ad_model,
      q, q_dirs,
      qdot, qdot_dirs,
      tau, tau_dirs,
      CS.QDDot_0, ad_CS.QDDot_0
    );
    // nominal code
    // ForwardDynamics(model, q, qdot, tau, CS.QDDot_0);
  }

  LOG << "=== Initial Loop Start ===" << std::endl;
  // we have to compute the standard accelerations first as we use them to
  // compute the effects of each test force
  for(ci = 0; ci < CS.size(); ci++) {

    {
      SUPPRESS_LOGGING;
      // derivative and nominal code
      UpdateKinematicsCustom(
        model, ad_model,
        NULL, NULL,
        NULL, NULL,
        &CS.QDDot_0, &ad_CS.QDDot_0
      );
      // nominal code
      // UpdateKinematicsCustom(model, NULL, NULL, &CS.QDDot_0);
    }

    if(CS.constraintType[ci] == ConstraintSet::ContactConstraint)
    {
      LOG << "body_id = " << CS.body[ci] << std::endl;
      LOG << "point = " << CS.point[ci] << std::endl;
      LOG << "normal = " << CS.normal[ci] << std::endl;
      LOG << "QDDot_0 = " << CS.QDDot_0.transpose() << std::endl;
      {
        SUPPRESS_LOGGING;
        // derivative and nominal code
//        ad_CS.point_accel_0[ci] = ad_CS.point_accel_0[ci].leftCols(ndirs).eval();
        CS.point_accel_0[ci] = AD::CalcPointAcceleration (
          model, ad_model,
          q, q_dirs,
          qdot, qdot_dirs,
          CS.QDDot_0, ad_CS.QDDot_0,
          CS.body[ci], CS.point[ci],
          ad_CS.point_accel_0[ci],
          false
        );
        for (unsigned int idir = 0; idir < ndirs; ++idir) {
          // NOTE CS.acceleration[ci] is a constant value
          ad_CS.a(ci, idir) = CS.normal[ci].dot(ad_CS.point_accel_0[ci].col(idir));
        }
        // nominal code
        // NOTE is already above
        // CS.point_accel_0[ci]
        //   = CalcPointAcceleration (
        //     model, Q, QDot , CS.QDDot_0, CS.body[ci], CS.point[ci], false
        //   );
        CS.a[ci] = - CS.acceleration[ci]
          + CS.normal[ci].dot(CS.point_accel_0[ci]);
      }
      LOG << "point_accel_0 = " << CS.point_accel_0[ci].transpose();
    }
    else
    {
      std::cerr << "Forward Dynamic Contact Kokkevis: unsupported constraint \
        type." << std::endl;
      assert(false);
      abort();
    }
  }

  // Now we can compute and apply the test forces and use their net effect
  // to compute the inverse articlated inertia to fill K.
  for (ci = 0; ci < CS.size(); ci++) {

    LOG << "=== Testforce Loop Start ===" << std::endl;

    unsigned int movable_body_id = 0;
    Vector3d point_global;
    std::vector<SpatialTransform> ad_X_tmp(ndirs, SpatialTransform());
    SpatialTransform X_tmp;
    SpatialVector v_tmp;

    switch (CS.constraintType[ci]) {

      case ConstraintSet::ContactConstraint:

        movable_body_id = GetMovableBodyId(model, CS.body[ci]);

        // assemble the test force
        LOG << "normal = " << CS.normal[ci].transpose() << std::endl;

        // derivative and nominal code
        point_global = AD::CalcBodyToBaseCoordinates(
          model, ad_model, q, q_dirs, CS.body[ci] , CS.point[ci],
          ad_CS.point_global, false
        );

        // nominal evaluation
        // NOTE evaluated in derivative code
        // point_global = CalcBodyToBaseCoordinates(
        //   model, Q, CS.body[ci], CS.point[ci], false
        // );

        LOG << "point_global = " << point_global.transpose() << std::endl;

        // derivative evaluation
        for (unsigned idir = 0; idir < ndirs; ++idir)
        {
          ad_X_tmp[idir].E = Matrix3d::Zero();
          ad_X_tmp[idir].r = -ad_CS.point_global.col(idir);
        }
        X_tmp.E = Matrix3d::Identity();
        X_tmp.r = -point_global;
        v_tmp = SpatialVector (
          0., 0., 0. , -CS.normal[ci][0], -CS.normal[ci][1], -CS.normal[ci][2]
        );
        applyAdjointSTSV (
          ndirs,
          X_tmp, ad_X_tmp,
          v_tmp,
          CS.f_t[ci],
          ad_CS.f_t[ci]
        );
        ad_CS.f_ext_constraints[movable_body_id] = ad_CS.f_t[ci];
        // nominal evaluation
        // CS.f_t[ci] = X_tmp.applyAdjoint(v_tmp);
        CS.f_ext_constraints[movable_body_id] = CS.f_t[ci];

        LOG << "f_t[" << movable_body_id << "] = " << CS.f_t[ci].transpose()
          << std::endl;

        {
          // derivative evaluation
          AD::ForwardDynamicsAccelerationDeltas(
            model, ad_model, CS, ad_CS,
            CS.QDDot_t, ad_CS.QDDot_t,
            movable_body_id,
            CS.f_ext_constraints,
            ad_CS.f_ext_constraints
          );
          // nominal evaluation
          // ForwardDynamicsAccelerationDeltas(
          //   model, CS, CS.QDDot_t , movable_body_id, CS.f_ext_constraints
          // );

          LOG << "QDDot_0 = " << CS.QDDot_0.transpose() << std::endl;
          LOG << "QDDot_t = " << (CS.QDDot_t + CS.QDDot_0).transpose()
            << std::endl;
          LOG << "QDDot_t - QDDot_0 = " << (CS.QDDot_t).transpose() << std::endl;
        }

        // derivative evaluation
        ad_CS.f_ext_constraints[movable_body_id].setZero();
        // nominal evaluation
        CS.f_ext_constraints[movable_body_id].setZero();

        // derivative evaluation
        ad_CS.QDDot_t += ad_CS.QDDot_0;
        // nominal evaluation
        CS.QDDot_t += CS.QDDot_0;

        // compute the resulting acceleration
        {
          SUPPRESS_LOGGING;
          // derivative evaluation
          AD::UpdateKinematicsCustom(
            model, ad_model,
            NULL, NULL,
            NULL, NULL,
            &CS.QDDot_t, &ad_CS.QDDot_t
          );
          // nominal evaluation
          // UpdateKinematicsCustom(model, NULL, NULL, &CS.QDDot_t);
        }

        // resize matrix once
//        ad_CS.point_accel_t = ad_CS.point_accel_t.leftCols(ndirs).eval();
        for(unsigned int cj = 0; cj < CS.size(); cj++)
        {
          { SUPPRESS_LOGGING;
            // derivative evaluation
            point_accel_t = AD::CalcPointAcceleration (
              model, ad_model,
              q, q_dirs,
              qdot, qdot_dirs,
              CS.QDDot_t, ad_CS.QDDot_t,
              CS.body[cj], CS.point[cj],
              ad_CS.point_accel_t,
              false
            );
            // nominal evaluation
            // point_accel_t = CalcPointAcceleration(
            //   model, q, qdot, CS.QDDot_t, CS.body[cj], CS.point[cj], false
            // );
          }

          LOG << "point_accel_0  = " << CS.point_accel_0[ci].transpose()
            << std::endl;
          LOG << "point_accel_t = " << point_accel_t.transpose() << std::endl;

          // derivative evaluation
          for (unsigned int idir = 0; idir < ndirs; ++idir) {
             ad_CS.K[idir](ci,cj)
              = CS.normal[cj].dot(
                ad_CS.point_accel_t.col(idir) - ad_CS.point_accel_0[cj].col(idir)
              );
           }
          // nominal evaluation
          CS.K(ci,cj) = CS.normal[cj].dot(point_accel_t - CS.point_accel_0[cj]);

        }

      break;

      default:

        std::cerr << "Forward Dynamic Contact Kokkevis: unsupported constraint \
          type." << std::endl;
        assert(false);
        abort();

      break;

    }

  }

  LOG << "K = " << std::endl << CS.K << std::endl;
  LOG << "a = " << std::endl << CS.a << std::endl;

#ifndef RBDL_USE_SIMPLE_MATH
  Eigen::HouseholderQR<MatrixNd> householderQR;
  Eigen::ColPivHouseholderQR<MatrixNd> colPivHouseholderQR;
  Eigen::PartialPivLU<MatrixNd> partialPivLU;
  switch (CS.linear_solver) {
    case (LinearSolverPartialPivLU) :
      // nominal evaluation
      partialPivLU.compute(CS.K);
      CS.force = partialPivLU.solve(CS.a);
      // derivative evaluation
      for (unsigned int idir = 0; idir < ndirs; ++idir) {
        ad_CS.force.col(idir) = partialPivLU.solve(
          ad_CS.a.col(idir) - ad_CS.K[idir]*CS.force
        );
      }
      break;
    case (LinearSolverColPivHouseholderQR) :
      // nominal evaluation
      colPivHouseholderQR.compute(CS.K);
      CS.force = colPivHouseholderQR.solve(CS.a);
      // derivative evaluation
      for (unsigned int idir = 0; idir < ndirs; ++idir) {
        ad_CS.force.col(idir) = colPivHouseholderQR.solve(
          ad_CS.a.col(idir) - ad_CS.K[idir]*CS.force
        );
      }
      break;
    case (LinearSolverHouseholderQR) :
      // nominal evaluation
      householderQR.compute(CS.K);
      CS.force = householderQR.solve(CS.a);
      // derivative evaluation
      for (unsigned int idir = 0; idir < ndirs; ++idir) {
        ad_CS.force.col(idir) = householderQR.solve(
          ad_CS.a.col(idir) - ad_CS.K[idir]*CS.force
        );
      }
      break;
    default:
      LOG << "Error: Invalid linear solver: " << CS.linear_solver << std::endl;
      assert (0);
      break;
  }
#else
  LOG << "Error: Linear solver not implemented for AD!" << std::endl;
  assert (0);
  bool solve_successful = LinSolveGaussElimPivot (CS.K, CS.a, CS.force);
  assert (solve_successful);
#endif



  LOG << "f = " << CS.force.transpose() << std::endl;

  for (ci = 0; ci < CS.size(); ci++) {
    unsigned int body_id = CS.body[ci];
    unsigned int movable_body_id = body_id;

    if (model.IsFixedBodyId(body_id)) {
      unsigned int fbody_id = body_id - model.fixed_body_discriminator;
      movable_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    }

    for (unsigned idir = 0; idir < ndirs; ++idir)
    {
      ad_CS.f_ext_constraints[movable_body_id].col(idir) -=
        ad_CS.f_t[ci].col(idir) * CS.force[ci]
          + CS.f_t[ci] * ad_CS.force(ci, idir);
    }

    CS.f_ext_constraints[movable_body_id] -= CS.f_t[ci] * CS.force[ci];
    LOG << "f_ext[" << movable_body_id << "] = " << CS.f_ext_constraints[movable_body_id].transpose() << std::endl;
  }

  {
    SUPPRESS_LOGGING;
    // derivative evaluation
    AD::ForwardDynamicsApplyConstraintForces (
      model, ad_model, tau, tau_dirs, CS, ad_CS, qddot, ad_qddot
    );
    // nominal evaluation
    // ForwardDynamicsApplyConstraintForces (model, tau, CS, qddot);
  }

  LOG << "QDDot after applying f_ext: " << qddot.transpose() << std::endl;
  return;

}

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
