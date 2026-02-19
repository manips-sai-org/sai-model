/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>

// #include "JointAD.h"
// #include "SpatialAlgebraOperatorsAD.h"

#include "DynamicsED.h"
#include "JointED.h"
#include "SpatialAlgebraOperatorsED.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

using namespace RigidBodyDynamics::Math;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace ED {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void InverseDynamics(
  RigidBodyDynamics::Model &model,
  EDModel &ed_model,
  Math::VectorNd const &q,
  Math::MatrixNd const &q_dirs,
  Math::VectorNd const &qdot,
  Math::MatrixNd const &qdot_dirs,
  Math::VectorNd const &qddot,
  Math::MatrixNd const &qddot_dirs,
  Math::VectorNd &tau,
  Math::MatrixNd &ed_tau,
  std::vector<Math::SpatialVector> const *f_ext,
  std::vector<std::vector<Math::SpatialVector> > const *f_ext_dirs
) {

  const unsigned int ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == qddot_dirs.cols());
  ed_model.resize_directions(ndirs);

  // Reset the velocity of the root body
  // nominal evaluation
  model.v[0].setZero ();
  model.a[0].set (0., 0., 0., -model.gravity[0], -model.gravity[1], -model.gravity[2]);
  // derivative evaluation
  ed_model.v[0].setZero ();
  ed_model.a[0].setZero();

  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    unsigned int q_index = model.mJoints[i].q_index;
    unsigned int lambda = model.lambda[i];

    jcalc (model, i, q, qdot);

    // nominal evaluation
    model.v[i] = model.X_lambda[i].apply(model.v[model.lambda[i]]);
    // derivative evaluation
    // d v[i] / d q
    ed_model.v[i].leftCols(ndirs)
        = crossm(model.v[i])*model.S[i]*q_dirs.row(model.mJoints[i].q_index)
        + model.X_lambda[i].toMatrix()*ed_model.v[lambda].leftCols(ndirs)
        + model.S[i]*qdot_dirs.row(model.mJoints[i].q_index);
    // nominal evaluation continued
    model.v[i] += model.v_J[i];

    // nominal evaluation
    model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);
    // derivative evaluation
    ed_model.c[i].leftCols(ndirs) =
        crossm(model.v[i])*model.S[i]*qdot_dirs.row(model.mJoints[i].q_index)
        - crossm(model.v_J[i]) * (ed_model.v[i].leftCols(ndirs) );

    if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom){
        // nominal evaluation
        model.a[i] = model.X_lambda[i].apply(model.a[model.lambda[i]]);
        // derivative evaluation
        ed_model.a[i] = crossm(model.a[i])*model.S[i]*q_dirs.row(model.mJoints[i].q_index)
            + model.X_lambda[i].toMatrix()*ed_model.a[lambda].leftCols(ndirs)
            + ed_model.c[i].leftCols(ndirs)
            + model.S[i] * qddot_dirs.row(model.mJoints[i].q_index).leftCols(ndirs);
        // nominal evaluation continued
        model.a[i] += model.c[i] + model.S[i] * qddot[q_index];
    } else if (model.mJoints[i].mDoFCount == 3) {
      cerr << "Multi-dof not supported." << endl;
      abort();
    } else if(model.mJoints[i].mJointType == JointTypeCustom){
      cerr << __FILE__ << " " << __LINE__
           << ": Custom joints not supported." << endl;
      abort();
    } else {
      cerr << __FILE__ << " " << __LINE__
           << ": Unknown unsupported joint." << endl;
      abort();
    }

    if (!model.mBodies[i].mIsVirtual) {
      // nominal evaluation
      ed_model.h[i] = model.I[i] * model.v[i];
      // derivative evaluation
      Math::MatrixNd const Ii_mat = model.I[i].toMatrix();

      // nominal evaluation
      model.f[i] = model.I[i] * model.a[i] + crossf(model.v[i], ed_model.h[i]);
      // derivative evaluation

      ed_model.f[i].leftCols(ndirs) =
          Ii_mat * (ed_model.a[i].leftCols(ndirs))
          + (crossf_rhs_T(ed_model.h[i])
          + crossf(model.v[i]) * Ii_mat) * (ed_model.v[i].leftCols(ndirs));
    } else {
      // nominal evaluation
      model.f[i].setZero();
      // derivative evaluation
      ed_model.f[i].setZero();
    }
  }

  if (f_ext != NULL) {
    cerr << __FILE__ << " " << __LINE__
         << ": External forces are not allowed." << endl;    abort();
    abort();
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {


    if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom){
        const unsigned int q_index = model.mJoints[i].q_index;
        // nominal evaluation
        tau[q_index] = model.S[i].dot(model.f[i]);
        // derivative evaluation
        ed_tau.row(q_index)
            = model.S[i].transpose() * ed_model.f[i].leftCols(ndirs);

    } else if (model.mJoints[i].mDoFCount == 3) {
      cerr << "Multi-dof not supported." << endl;
      abort();
    } else if (model.mJoints[i].mJointType == JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__
           << ": Custom joints not supported." << endl;    abort();
      abort();
    } else {
      cerr << __FILE__ << " " << __LINE__
           << ": Unknown unsupported joint." << endl;
      abort();
    }

    if (model.lambda[i] != 0) {
      // nominal evaluation
      model.f[model.lambda[i]] += model.X_lambda[i].applyTranspose(model.f[i]);
      // derivative evaluation
      // d a[i] / d q
      ed_model.f[model.lambda[i]].leftCols(ndirs) +=
          model.X_lambda[i].toMatrixTranspose()
          * (crossf_rhs_T(model.f[i]) * model.S[i]*q_dirs.row(model.mJoints[i].q_index).leftCols(ndirs)
          + ed_model.f[i].leftCols(ndirs));
    }
  }

  return;
}

RBDL_DLLAPI void NonlinearEffects (
    Model & model,
    EDModel & ed_model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot,
    const Math::MatrixNd & qdot_dirs,
    Math::VectorNd & tau,
    Math::MatrixNd & ed_tau
) {
  // LOG << "-------- " << __func__ << " --------" << std::endl;
  const unsigned int ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  // assert(ndirs == qddot_dirs.cols());
  ed_model.resize_directions(ndirs);

  SpatialVector spatial_gravity (0., 0., 0., -model.gravity[0], -model.gravity[1], -model.gravity[2]);

  // Reset the velocity of the root body
  // nominal evaluation
  model.v[0].setZero();
  model.a[0] = spatial_gravity;

  // derivative evaluation
  ed_model.v[0].leftCols(ndirs).setZero();
  ed_model.a[0].leftCols(ndirs).setZero();

  for (unsigned int i = 1; i < model.mJointUpdateOrder.size(); i++) {
    jcalc (model, model.mJointUpdateOrder[i], q, qdot);
  }

  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    if (model.lambda[i] == 0) {
      // nominal evaluation
      model.v[i] = model.v_J[i];
      // derivative evaluation
      // NOTE bodyidx2s1idx[i] gives index of 1 entry in model.S[i]
      ed_model.v[i].row(ed_model.bodyidx2s1idx[i]).leftCols(ndirs)
        = qdot_dirs.row(model.mJoints[i].q_index);

      // nominal evaluation
      model.a[i] = model.X_lambda[i].apply(spatial_gravity);
      // derivative evaluation
      ed_model.a[i].leftCols(ndirs)
        = crossm(model.a[i]).col(ed_model.bodyidx2s1idx[i])
        * q_dirs.row(model.mJoints[i].q_index);
    } else {
      // nominal evaluation
      model.v[i] = model.X_lambda[i].apply(model.v[model.lambda[i]]);
      // derivative evaluation
      // NOTE we compute in 4 steps to not evaluate zero blocks of X.apply(dirs)
      ed_model.v[i].leftCols(ndirs)
          = crossm(model.v[i]).col(ed_model.bodyidx2s1idx[i])
          * q_dirs.row(model.mJoints[i].q_index);
      ed_model.v[i].block(0, 0, 3, ndirs)
        += model.X_lambda[i].E*ed_model.v[model.lambda[i]].block(0, 0, 3, ndirs);
      ed_model.v[i].block(3, 0, 3, ndirs) -= model.X_lambda[i].E * (
          VectorCrossMatrix(model.X_lambda[i].r)
          * ed_model.v[model.lambda[i]].block(0, 0, 3, ndirs)
          - ed_model.v[model.lambda[i]].block(3, 0, 3, ndirs)
        );
      ed_model.v[i].row(ed_model.bodyidx2s1idx[i]).leftCols(ndirs)
          += qdot_dirs.row(model.mJoints[i].q_index);
      // nominal evaluation continued
      model.v[i] += model.v_J[i];

      // nominal evaluation
      model.c[i] = model.c_J[i] + crossm(model.v[i], model.v_J[i]);
      // derivative evaluation
      // NOTE we compute in 3 steps to save zero block evaluation
      ed_model.c[i].leftCols(ndirs)
           = crossm(model.v[i]).col(ed_model.bodyidx2s1idx[i])*qdot_dirs.row(model.mJoints[i].q_index);
      ed_model.c[i].block(0, 0, 3, ndirs)
          -= VectorCrossMatrix(model.v_J[i].head<3>())
          * ed_model.v[i].block(0, 0, 3, ndirs);
      ed_model.c[i].block(3, 0, 3, ndirs)
          -= VectorCrossMatrix(model.v_J[i].tail<3>()) * (ed_model.v[i].block(0, 0, 3, ndirs))
          + VectorCrossMatrix(model.v_J[i].head<3>()) * (ed_model.v[i].block(3, 0, 3, ndirs));

      // nominal evaluation
      model.a[i] = model.X_lambda[i].apply(model.a[model.lambda[i]]);
      // derivative evaluation
      ed_model.a[i].leftCols(ndirs) =
          crossm(model.a[i]).col(ed_model.bodyidx2s1idx[i])*q_dirs.row(model.mJoints[i].q_index)
          + ed_model.c[i].leftCols(ndirs);
      ed_model.a[i].block(0, 0, 3, ndirs)
        += model.X_lambda[i].E*ed_model.a[model.lambda[i]].block(0, 0, 3, ndirs);
      ed_model.a[i].block(3, 0, 3, ndirs)
        -= model.X_lambda[i].E * (
          VectorCrossMatrix(model.X_lambda[i].r)
          * ed_model.a[model.lambda[i]].block(0, 0, 3, ndirs)
          - ed_model.a[model.lambda[i]].block(3, 0, 3, ndirs)
        );
      // nominal evaluation continued
      model.a[i] += model.c[i];
    }

    if (!model.mBodies[i].mIsVirtual) {
      // nominal evaluation
      ed_model.h[i] = model.I[i] * model.v[i];
      Math::SpatialMatrix const Ii_mat = model.I[i].toMatrix();

      SpatialMatrix cross_rhs_hi_T = crossf_rhs_T(ed_model.h[i]);
      SpatialMatrix cross_lhs_vi = crossf(model.v[i]);

      // nominal evaluation
      model.f[i] = model.I[i] * model.a[i] + crossf(model.v[i], ed_model.h[i]);
      // derivative evaluation
      // TODO save zero block evaluation here
      ed_model.f[i].leftCols(ndirs)
          = Ii_mat*(ed_model.a[i].leftCols(ndirs))
          + (cross_rhs_hi_T + cross_lhs_vi * Ii_mat) * (ed_model.v[i].leftCols(ndirs));
    } else {
      model.f[i].setZero();
      ed_model.f[i].leftCols(ndirs).setZero();
    }
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
    if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom){
        const unsigned int q_index = model.mJoints[i].q_index;
        // nominal evaluation
        tau[q_index] = model.f[i](ed_model.bodyidx2s1idx[i]);
      // derivative evaluation
      // d tau [i] = d tau [i] / d q + d tau [i] / d qdot
        ed_tau.leftCols(ndirs).row(q_index)
          = ed_model.f[i].leftCols(ndirs).row(ed_model.bodyidx2s1idx[i]);

    } else if(model.mJoints[i].mJointType == JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__
           << ": Custom joints not supported." << endl;    abort();
      abort();
    } else if (model.mJoints[i].mDoFCount == 3) {
      cerr << "Multi-dof not supported." << endl;
      abort();
    } else {
      cerr << __FILE__ << " " << __LINE__
           << ": Unknown unsupported joint." << endl;
      abort();
    }

    if (model.lambda[i] != 0) {
      // nominal evaluation
      model.f[model.lambda[i]] += model.X_lambda[i].applyTranspose(model.f[i]);
      // derivative evaluation
      // TODO save zero block computation here
      ed_model.f[model.lambda[i]].leftCols(ndirs)
          += model.X_lambda[i].toMatrixTranspose()
          * (crossf_rhs_T(model.f[i]).col(ed_model.bodyidx2s1idx[i])*q_dirs.row(model.mJoints[i].q_index).leftCols(ndirs)
             + ed_model.f[i].leftCols(ndirs));
    }
  }
}

RBDL_DLLAPI
void CompositeRigidBodyAlgorithm (
  Model &model,
  EDModel &ed_model,
  Math::VectorNd const & q,
  Math::MatrixNd const & q_dirs,
  Math::MatrixNd & H,
  std::vector<Math::MatrixNd> & H_dirs,
  bool update_kinematics
)
{
  assert (H.rows() == model.dof_count && H.cols() == model.dof_count);

  // get number if directions
  const size_t ndirs = q_dirs.cols();
  ed_model.resize_directions(ndirs);

  if (update_kinematics)
  {
    for (unsigned int i = 1; i < model.mBodies.size(); i++)
    {
      jcalc_X_lambda_S (model, i, q);
    }
  }

  for (unsigned int i = 1; i < model.mBodies.size(); i++)
  {
    // nominal evaluation
    model.Ic[i] = model.I[i];
    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; idir++)
    {
      ed_model.Ic[i][idir].setZero();
    }
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--)
  {
    if (model.lambda[i] != 0)
    {
      // nominal evaluation
      // NOTE we require temporary spatial rbi for efficient computation
      SpatialRigidBodyInertia temp = model.X_lambda[i].applyTranspose(model.Ic[i]);
      model.Ic[model.lambda[i]] += temp;

      // derivative evaluation
      // NOTE we have to transform back S vector
      const Math::SpatialVector imv = model.X_lambda[i].inverse().apply(model.S[i]);

      for (unsigned idir = 0; idir < ndirs; idir++) {
        Vector3d E_T_mr = model.X_lambda[i].E.transpose() * ed_model.Ic[i][idir].h;
        SpatialRigidBodyInertia rbi = SpatialRigidBodyInertia (
          0.,
          E_T_mr,
          model.X_lambda[i].E.transpose() *
          Matrix3d (
            ed_model.Ic[i][idir].Ixx, ed_model.Ic[i][idir].Iyx, ed_model.Ic[i][idir].Izx,
            ed_model.Ic[i][idir].Iyx, ed_model.Ic[i][idir].Iyy, ed_model.Ic[i][idir].Izy,
            ed_model.Ic[i][idir].Izx, ed_model.Ic[i][idir].Izy, ed_model.Ic[i][idir].Izz
            ) * model.X_lambda[i].E
          - VectorCrossMatrix(model.X_lambda[i].r) * VectorCrossMatrix (model.X_lambda[i].E.transpose() * ed_model.Ic[i][idir].h)
          - VectorCrossMatrix (E_T_mr) * VectorCrossMatrix (model.X_lambda[i].r)
        );

        const Math::Vector3d w = imv.head(3)*q_dirs(model.mJoints[i].q_index, idir);
        const Math::Vector3d v0 = imv.tail(3)*q_dirs(model.mJoints[i].q_index, idir);
        ed_model.Ic[model.lambda[i]][idir]
          += SpatialRigidBodyInertia (
              0,
              w.cross(temp.h) + temp.m * v0,
              Matrix3d(
                -temp.Iyx*w[2] + temp.Izx*w[1] - temp.Iyx*w[2] + temp.Izx*w[1] + 2.*(temp.h[1]*v0[1] + temp.h[2]*v0[2]),
                 temp.Ixx*w[2] - temp.Izx*w[0] - temp.Iyy*w[2] + temp.Izy*w[1] -     temp.h[0]*v0[1] - temp.h[1]*v0[0] ,
                -temp.Ixx*w[1] + temp.Iyx*w[0] - temp.Izy*w[2] + temp.Izz*w[1] -     temp.h[0]*v0[2] - temp.h[2]*v0[0] ,

                 temp.Ixx*w[2] - temp.Iyy*w[2] + temp.Izy*w[1] - temp.Izx*w[0] -     temp.h[0]*v0[1] - temp.h[1]*v0[0] ,
                 temp.Iyx*w[2] + temp.Iyx*w[2] - temp.Izy*w[0] - temp.Izy*w[0] + 2.*(temp.h[0]*v0[0] + temp.h[2]*v0[2]),
                 temp.Izx*w[2] - temp.Iyx*w[1] + temp.Iyy*w[0] - temp.Izz*w[0] -     temp.h[1]*v0[2] - temp.h[2]*v0[1] ,

                -temp.Ixx*w[1] + temp.Iyx*w[0] - temp.Izy*w[2] + temp.Izz*w[1] -     temp.h[0]*v0[2] - temp.h[2]*v0[0] ,
                -temp.Iyx*w[1] + temp.Iyy*w[0] + temp.Izx*w[2] - temp.Izz*w[0] -     temp.h[1]*v0[2] - temp.h[2]*v0[1] ,
                -temp.Izx*w[1] + temp.Izy*w[0] - temp.Izx*w[1] + temp.Izy*w[0] + 2.*(temp.h[0]*v0[0] + temp.h[1]*v0[1])
              )
          )
          + rbi
        ;
      }
    }

    unsigned int dof_index_i = model.mJoints[i].q_index;

    if (model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom)
    {
      // nominal evaluation
      SpatialVector F             = model.Ic[i] * model.S[i];
      H(dof_index_i, dof_index_i) = model.S[i].dot(F);

      // derivative evaluation
      for (unsigned idir = 0; idir < ndirs; idir++)
      {
        ed_model.F[i].col(idir) = ed_model.Ic[i][idir] * model.S[i];
        H_dirs[idir](dof_index_i, dof_index_i) = model.S[i].dot(ed_model.F[i].col(idir));
      }

      unsigned int j = i;
      unsigned int dof_index_j = dof_index_i;

      while (model.lambda[j] != 0) {
        // derivative evaluation
        // TODO do not evaluate zero blocks
        ed_model.F[i].leftCols(ndirs) = model.X_lambda[j].toMatrixTranspose()
            * (ed_model.F[i].leftCols(ndirs) + crossf_rhs_T(F)*model.S[j]*q_dirs.row(model.mJoints[j].q_index).leftCols(ndirs));
        // nominal evaluation
        F = model.X_lambda[j].applyTranspose(F);

        j = model.lambda[j];
        dof_index_j = model.mJoints[j].q_index;

        // nominal evaluation
        H(dof_index_i,dof_index_j) = F.dot(model.S[j]);
        H(dof_index_j,dof_index_i) = H(dof_index_i,dof_index_j);
        // derivative evaluation
        for (unsigned idir = 0; idir < ndirs; idir++) {
          H_dirs[idir](dof_index_i, dof_index_j)
            = ed_model.F[i].col(idir).dot(model.S[j]);
          H_dirs[idir](dof_index_j, dof_index_i)
            = H_dirs[idir](dof_index_i, dof_index_j);
        }
      }
    } else if (model.mJoints[i].mDoFCount == 3 && model.mJoints[i].mJointType != JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__ << ":"
           << "Multi-DoF joint not supported." << endl;
      abort();
    } else if (model.mJoints[i].mJointType == JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__ << ":"
           << " Custom joints not supported." << endl;
      abort();
    } else {
      cerr << __FILE__ << " " << __LINE__
           << ": Unknown unsupported joint." << endl;
      abort();
    }
  }
}

// -----------------------------------------------------------------------------
} // ED
// -----------------------------------------------------------------------------
} // RigidBodyDynamics
// -----------------------------------------------------------------------------
