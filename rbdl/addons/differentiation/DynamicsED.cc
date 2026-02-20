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

// // original
// RBDL_DLLAPI void NonlinearEffects (
//     Model & model,
//     EDModel & ed_model,
//     const Math::VectorNd & q,
//     const Math::MatrixNd & q_dirs,
//     const Math::VectorNd & qdot,
//     const Math::MatrixNd & qdot_dirs,
//     Math::VectorNd & tau,
//     Math::MatrixNd & ed_tau
// ) {
//   // LOG << "-------- " << __func__ << " --------" << std::endl;
//   const unsigned int ndirs = q_dirs.cols();
//   assert(ndirs == qdot_dirs.cols());
//   // assert(ndirs == qddot_dirs.cols());
//   ed_model.resize_directions(ndirs);

//   SpatialVector spatial_gravity (0., 0., 0., -model.gravity[0], -model.gravity[1], -model.gravity[2]);

//   // Reset the velocity of the root body
//   // nominal evaluation
//   model.v[0].setZero();
//   model.a[0] = spatial_gravity;

//   // derivative evaluation
//   ed_model.v[0].leftCols(ndirs).setZero();
//   ed_model.a[0].leftCols(ndirs).setZero();

//   for (unsigned int i = 1; i < model.mJointUpdateOrder.size(); i++) {
//     jcalc (model, model.mJointUpdateOrder[i], q, qdot);
//   }

//   for (unsigned int i = 1; i < model.mBodies.size(); i++) {
//     if (model.lambda[i] == 0) {
//       // nominal evaluation
//       model.v[i] = model.v_J[i];
//       // derivative evaluation
//       // NOTE bodyidx2s1idx[i] gives index of 1 entry in model.S[i]
//       ed_model.v[i].row(ed_model.bodyidx2s1idx[i]).leftCols(ndirs)
//         = qdot_dirs.row(model.mJoints[i].q_index);

//       // nominal evaluation
//       model.a[i] = model.X_lambda[i].apply(spatial_gravity);
//       // derivative evaluation
//       ed_model.a[i].leftCols(ndirs)
//         = crossm(model.a[i]).col(ed_model.bodyidx2s1idx[i])
//         * q_dirs.row(model.mJoints[i].q_index);
//     } else {
//       // nominal evaluation
//       model.v[i] = model.X_lambda[i].apply(model.v[model.lambda[i]]);
//       // derivative evaluation
//       // NOTE we compute in 4 steps to not evaluate zero blocks of X.apply(dirs)
//       ed_model.v[i].leftCols(ndirs)
//           = crossm(model.v[i]).col(ed_model.bodyidx2s1idx[i])
//           * q_dirs.row(model.mJoints[i].q_index);
//       ed_model.v[i].block(0, 0, 3, ndirs)
//         += model.X_lambda[i].E*ed_model.v[model.lambda[i]].block(0, 0, 3, ndirs);
//       ed_model.v[i].block(3, 0, 3, ndirs) -= model.X_lambda[i].E * (
//           VectorCrossMatrix(model.X_lambda[i].r)
//           * ed_model.v[model.lambda[i]].block(0, 0, 3, ndirs)
//           - ed_model.v[model.lambda[i]].block(3, 0, 3, ndirs)
//         );
//       ed_model.v[i].row(ed_model.bodyidx2s1idx[i]).leftCols(ndirs)
//           += qdot_dirs.row(model.mJoints[i].q_index);
//       // nominal evaluation continued
//       model.v[i] += model.v_J[i];

//       // nominal evaluation
//       model.c[i] = model.c_J[i] + crossm(model.v[i], model.v_J[i]);
//       // derivative evaluation
//       // NOTE we compute in 3 steps to save zero block evaluation
//       ed_model.c[i].leftCols(ndirs)
//            = crossm(model.v[i]).col(ed_model.bodyidx2s1idx[i])*qdot_dirs.row(model.mJoints[i].q_index);
//       ed_model.c[i].block(0, 0, 3, ndirs)
//           -= VectorCrossMatrix(model.v_J[i].head<3>())
//           * ed_model.v[i].block(0, 0, 3, ndirs);
//       ed_model.c[i].block(3, 0, 3, ndirs)
//           -= VectorCrossMatrix(model.v_J[i].tail<3>()) * (ed_model.v[i].block(0, 0, 3, ndirs))
//           + VectorCrossMatrix(model.v_J[i].head<3>()) * (ed_model.v[i].block(3, 0, 3, ndirs));

//       // nominal evaluation
//       model.a[i] = model.X_lambda[i].apply(model.a[model.lambda[i]]);
//       // derivative evaluation
//       ed_model.a[i].leftCols(ndirs) =
//           crossm(model.a[i]).col(ed_model.bodyidx2s1idx[i])*q_dirs.row(model.mJoints[i].q_index)
//           + ed_model.c[i].leftCols(ndirs);
//       ed_model.a[i].block(0, 0, 3, ndirs)
//         += model.X_lambda[i].E*ed_model.a[model.lambda[i]].block(0, 0, 3, ndirs);
//       ed_model.a[i].block(3, 0, 3, ndirs)
//         -= model.X_lambda[i].E * (
//           VectorCrossMatrix(model.X_lambda[i].r)
//           * ed_model.a[model.lambda[i]].block(0, 0, 3, ndirs)
//           - ed_model.a[model.lambda[i]].block(3, 0, 3, ndirs)
//         );
//       // nominal evaluation continued
//       model.a[i] += model.c[i];
//     }

//     if (!model.mBodies[i].mIsVirtual) {
//       // nominal evaluation
//       ed_model.h[i] = model.I[i] * model.v[i];
//       Math::SpatialMatrix const Ii_mat = model.I[i].toMatrix();

//       SpatialMatrix cross_rhs_hi_T = crossf_rhs_T(ed_model.h[i]);
//       SpatialMatrix cross_lhs_vi = crossf(model.v[i]);

//       // nominal evaluation
//       model.f[i] = model.I[i] * model.a[i] + crossf(model.v[i], ed_model.h[i]);
//       // derivative evaluation
//       // TODO save zero block evaluation here
//       ed_model.f[i].leftCols(ndirs)
//           = Ii_mat*(ed_model.a[i].leftCols(ndirs))
//           + (cross_rhs_hi_T + cross_lhs_vi * Ii_mat) * (ed_model.v[i].leftCols(ndirs));
//     } else {
//       model.f[i].setZero();
//       ed_model.f[i].leftCols(ndirs).setZero();
//     }
//   }

//   for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
//     if(model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom){
//         const unsigned int q_index = model.mJoints[i].q_index;
//         // nominal evaluation
//         tau[q_index] = model.f[i](ed_model.bodyidx2s1idx[i]);
//       // derivative evaluation
//       // d tau [i] = d tau [i] / d q + d tau [i] / d qdot
//         ed_tau.leftCols(ndirs).row(q_index)
//           = ed_model.f[i].leftCols(ndirs).row(ed_model.bodyidx2s1idx[i]);

//     } else if(model.mJoints[i].mJointType == JointTypeCustom) {
//       cerr << __FILE__ << " " << __LINE__
//            << ": Custom joints not supported." << endl;    abort();
//       abort();
//     } else if (model.mJoints[i].mDoFCount == 3) {
//       cerr << "Multi-dof not supported." << endl;
//       abort();
//     } else {
//       cerr << __FILE__ << " " << __LINE__
//            << ": Unknown unsupported joint." << endl;
//       abort();
//     }

//     if (model.lambda[i] != 0) {
//       // nominal evaluation
//       model.f[model.lambda[i]] += model.X_lambda[i].applyTranspose(model.f[i]);
//       // derivative evaluation
//       // TODO save zero block computation here
//       ed_model.f[model.lambda[i]].leftCols(ndirs)
//           += model.X_lambda[i].toMatrixTranspose()
//           * (crossf_rhs_T(model.f[i]).col(ed_model.bodyidx2s1idx[i])*q_dirs.row(model.mJoints[i].q_index).leftCols(ndirs)
//              + ed_model.f[i].leftCols(ndirs));
//     }
//   }
// }

// // optimized
// RBDL_DLLAPI void NonlinearEffects (
//     Model & model,
//     EDModel & ed_model,
//     const Math::VectorNd & q,
//     const Math::MatrixNd & q_dirs,
//     const Math::VectorNd & qdot,
//     const Math::MatrixNd & qdot_dirs,
//     Math::VectorNd & tau,
//     Math::MatrixNd & ed_tau
// ) {
//   const unsigned int ndirs = q_dirs.cols();
//   const unsigned int num_bodies = model.mBodies.size();
//   ed_model.resize_directions(ndirs);

//   // Pre-define spatial gravity to avoid re-construction
//   const SpatialVector spatial_gravity (0., 0., 0., -model.gravity[0], -model.gravity[1], -model.gravity[2]);

//   model.v[0].setZero();
//   model.a[0] = spatial_gravity;
//   ed_model.v[0].leftCols(ndirs).setZero();
//   ed_model.a[0].leftCols(ndirs).setZero();

//   // Forward Pass: Kinematics and Derivatives
//   for (unsigned int i = 1; i < num_bodies; i++) {
//     jcalc (model, i, q, qdot); // Direct indexing is faster than mJointUpdateOrder if tree is linear
    
//     const unsigned int lambda = model.lambda[i];
//     const unsigned int q_idx = model.mJoints[i].q_index;
//     const unsigned int s1_idx = ed_model.bodyidx2s1idx[i];
//     SpatialTransform X_lam = model.X_lambda[i];

//     if (lambda == 0) {
//       model.v[i] = model.v_J[i];
//       ed_model.v[i].leftCols(ndirs).setZero();
//       ed_model.v[i].row(s1_idx).leftCols(ndirs) = qdot_dirs.row(q_idx);

//       model.a[i] = X_lam.apply(spatial_gravity);
//       // Pre-extract column of cross-matrix for q_dirs multiplication
//       ed_model.a[i].leftCols(ndirs) = crossm(model.a[i]).col(s1_idx) * q_dirs.row(q_idx);
//     } else {
//       // Velocity nominal and derivative
//       model.v[i] = X_lam.apply(model.v[lambda]);
      
//       // Optimization: Group block operations to improve cache locality
//       auto v_curr = ed_model.v[i].leftCols(ndirs);
//       auto v_lam  = ed_model.v[lambda].leftCols(ndirs);
      
//       v_curr = crossm(model.v[i]).col(s1_idx) * q_dirs.row(q_idx);
//       v_curr.block(0, 0, 3, ndirs) += X_lam.E * v_lam.block(0, 0, 3, ndirs);
//       v_curr.block(3, 0, 3, ndirs) -= X_lam.E * (VectorCrossMatrix(X_lam.r) * v_lam.block(0, 0, 3, ndirs) - v_lam.block(3, 0, 3, ndirs));
//       v_curr.row(s1_idx) += qdot_dirs.row(q_idx);
      
//       model.v[i] += model.v_J[i];

//       // Acceleration center (c_i) derivative
//       model.c[i] = model.c_J[i] + crossm(model.v[i], model.v_J[i]);
      
//       auto c_curr = ed_model.c[i].leftCols(ndirs);
//       const Vector3d& vJ_head = model.v_J[i].head<3>();
//       const Vector3d& vJ_tail = model.v_J[i].tail<3>();

//       c_curr = crossm(model.v[i]).col(s1_idx) * qdot_dirs.row(q_idx);
//       c_curr.block(0, 0, 3, ndirs) -= VectorCrossMatrix(vJ_head) * v_curr.block(0, 0, 3, ndirs);
//       c_curr.block(3, 0, 3, ndirs) -= (VectorCrossMatrix(vJ_tail) * v_curr.block(0, 0, 3, ndirs) + VectorCrossMatrix(vJ_head) * v_curr.block(3, 0, 3, ndirs));

//       // Acceleration nominal and derivative
//       model.a[i] = X_lam.apply(model.a[lambda]);
      
//       auto a_curr = ed_model.a[i].leftCols(ndirs);
//       auto a_lam  = ed_model.a[lambda].leftCols(ndirs);
      
//       a_curr = crossm(model.a[i]).col(s1_idx) * q_dirs.row(q_idx) + c_curr;
//       a_curr.block(0, 0, 3, ndirs) += X_lam.E * a_lam.block(0, 0, 3, ndirs);
//       a_curr.block(3, 0, 3, ndirs) -= X_lam.E * (VectorCrossMatrix(X_lam.r) * a_lam.block(0, 0, 3, ndirs) - a_lam.block(3, 0, 3, ndirs));
      
//       model.a[i] += model.c[i];
//     }

//     // Force evaluation
//     if (!model.mBodies[i].mIsVirtual) {
//       ed_model.h[i] = model.I[i] * model.v[i];
//       const SpatialMatrix Ii_mat = model.I[i].toMatrix();
      
//       // Cache crossf terms
//       const SpatialMatrix cross_combined = crossf_rhs_T(ed_model.h[i]) + crossf(model.v[i]) * Ii_mat;
      
//       model.f[i] = model.I[i] * model.a[i] + crossf(model.v[i], ed_model.h[i]);
//       ed_model.f[i].leftCols(ndirs) = Ii_mat * ed_model.a[i].leftCols(ndirs) + cross_combined * ed_model.v[i].leftCols(ndirs);
//     } else {
//       model.f[i].setZero();
//       ed_model.f[i].leftCols(ndirs).setZero();
//     }
//   }

//   // Backward Pass: Torques and Force Propagation
//   for (unsigned int i = num_bodies - 1; i > 0; i--) {
//     const unsigned int q_idx = model.mJoints[i].q_index;
//     const unsigned int s1_idx = ed_model.bodyidx2s1idx[i];
//     const unsigned int lambda = model.lambda[i];

//     // Joint support check
//     if (model.mJoints[i].mDoFCount != 1 || model.mJoints[i].mJointType == JointTypeCustom) {
//       std::cerr << "Unsupported joint type at body " << i << std::endl;
//       abort();
//     }

//     tau[q_idx] = model.f[i](s1_idx);
//     ed_tau.row(q_idx).leftCols(ndirs) = ed_model.f[i].row(s1_idx).leftCols(ndirs);

//     if (lambda != 0) {
//       model.f[lambda] += model.X_lambda[i].applyTranspose(model.f[i]);
      
//       // Optimization: Compute crossf_rhs_T * q_dirs once then add to ed_model.f
//       ed_model.f[lambda].leftCols(ndirs) += model.X_lambda[i].toMatrixTranspose() * (
//           crossf_rhs_T(model.f[i]).col(s1_idx) * q_dirs.row(q_idx) + ed_model.f[i].leftCols(ndirs)
//       );
//     }
//   }
// }

// only derivatives
RBDL_DLLAPI void NonlinearEffects (
    Model & model,
    EDModel & ed_model,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs,
    const Math::VectorNd & qdot,
    const Math::MatrixNd & qdot_dirs,
    Math::MatrixNd & ed_tau,
    const bool gravity_opt
) {
  const unsigned int ndirs = q_dirs.cols();
  const unsigned int num_bodies = model.mBodies.size();
  ed_model.resize_directions(ndirs);

  SpatialVector spatial_gravity (0., 0., 0., -model.gravity[0], -model.gravity[1], -model.gravity[2]);
  if (!gravity_opt) {
    spatial_gravity.setZero();
  }

  model.v[0].setZero();
  model.a[0] = spatial_gravity;
  ed_model.v[0].leftCols(ndirs).setZero();
  ed_model.a[0].leftCols(ndirs).setZero();

  // Forward Pass
  for (unsigned int i = 1; i < num_bodies; i++) {
    jcalc (model, i, q, qdot); 
    
    const unsigned int lambda = model.lambda[i];
    const unsigned int q_idx = model.mJoints[i].q_index;
    const unsigned int s1_idx = ed_model.bodyidx2s1idx[i];
    SpatialTransform X_lam = model.X_lambda[i];

    if (lambda == 0) {
      model.v[i] = model.v_J[i];
      ed_model.v[i].leftCols(ndirs).setZero();
      ed_model.v[i].row(s1_idx).leftCols(ndirs) = qdot_dirs.row(q_idx);

      model.a[i] = X_lam.apply(spatial_gravity);
      ed_model.a[i].leftCols(ndirs) = crossm(model.a[i]).col(s1_idx) * q_dirs.row(q_idx);
    } else {
      model.v[i] = X_lam.apply(model.v[lambda]);
      
      auto v_curr = ed_model.v[i].leftCols(ndirs);
      auto v_lam  = ed_model.v[lambda].leftCols(ndirs);
      
      v_curr = crossm(model.v[i]).col(s1_idx) * q_dirs.row(q_idx);
      v_curr.block(0, 0, 3, ndirs) += X_lam.E * v_lam.block(0, 0, 3, ndirs);
      v_curr.block(3, 0, 3, ndirs) -= X_lam.E * (VectorCrossMatrix(X_lam.r) * v_lam.block(0, 0, 3, ndirs) - v_lam.block(3, 0, 3, ndirs));
      v_curr.row(s1_idx) += qdot_dirs.row(q_idx);
      
      model.v[i] += model.v_J[i];

      model.c[i] = model.c_J[i] + crossm(model.v[i], model.v_J[i]);
      
      auto c_curr = ed_model.c[i].leftCols(ndirs);
      const Vector3d& vJ_head = model.v_J[i].head<3>();
      const Vector3d& vJ_tail = model.v_J[i].tail<3>();

      c_curr = crossm(model.v[i]).col(s1_idx) * qdot_dirs.row(q_idx);
      c_curr.block(0, 0, 3, ndirs) -= VectorCrossMatrix(vJ_head) * v_curr.block(0, 0, 3, ndirs);
      c_curr.block(3, 0, 3, ndirs) -= (VectorCrossMatrix(vJ_tail) * v_curr.block(0, 0, 3, ndirs) + VectorCrossMatrix(vJ_head) * v_curr.block(3, 0, 3, ndirs));

      model.a[i] = X_lam.apply(model.a[lambda]);
      
      auto a_curr = ed_model.a[i].leftCols(ndirs);
      auto a_lam  = ed_model.a[lambda].leftCols(ndirs);
      
      a_curr = crossm(model.a[i]).col(s1_idx) * q_dirs.row(q_idx) + c_curr;
      a_curr.block(0, 0, 3, ndirs) += X_lam.E * a_lam.block(0, 0, 3, ndirs);
      a_curr.block(3, 0, 3, ndirs) -= X_lam.E * (VectorCrossMatrix(X_lam.r) * a_lam.block(0, 0, 3, ndirs) - a_lam.block(3, 0, 3, ndirs));
      
      model.a[i] += model.c[i];
    }

    if (!model.mBodies[i].mIsVirtual) {
      SpatialRigidBodyInertia Ii = model.I[i];
      ed_model.h[i] = Ii * model.v[i];
      model.f[i] = Ii * model.a[i] + crossf(model.v[i], ed_model.h[i]);

      // OPTIMIZATION: Replacing crossf(model.v[i]) * Ii_mat * ed_v with crossf(v, I*dv)
      // and crossf_rhs_T(h) * ed_v with crossf(dv, h)
      for (unsigned idir = 0; idir < ndirs; idir++) {
        SpatialVector dv = ed_model.v[i].col(idir);
        ed_model.f[i].col(idir) = Ii * ed_model.a[i].col(idir) 
                                + crossf(dv, ed_model.h[i]) 
                                + crossf(model.v[i], Ii * dv);
      }
    } else {
      model.f[i].setZero();
      ed_model.f[i].leftCols(ndirs).setZero();
    }
  }

  // Backward Pass
  for (unsigned int i = num_bodies - 1; i > 0; i--) {
    const unsigned int q_idx = model.mJoints[i].q_index;
    const unsigned int s1_idx = ed_model.bodyidx2s1idx[i];
    const unsigned int lambda = model.lambda[i];

    if (model.mJoints[i].mDoFCount != 1 || model.mJoints[i].mJointType == JointTypeCustom) {
      std::cerr << "Unsupported joint type at body " << i << std::endl;
      abort();
    }

    // tau removed.
    ed_tau.row(q_idx).leftCols(ndirs) = ed_model.f[i].row(s1_idx).leftCols(ndirs);

    if (lambda != 0) {
      model.f[lambda] += model.X_lambda[i].applyTranspose(model.f[i]);
      
      // OPTIMIZATION: Vectorized cross product and applyTranspose instead of 6x6 MatrixTranspose
      // Equivalent to: X^T * (ed_f + (S x* f) * q_dir)
      SpatialVector f_cross_S = crossf(model.S[i], model.f[i]);

      for (unsigned idir = 0; idir < ndirs; idir++) {
          SpatialVector ed_f_i_combined = ed_model.f[i].col(idir) 
                                        + f_cross_S * q_dirs(q_idx, idir);
          
          ed_model.f[lambda].col(idir) += model.X_lambda[i].applyTranspose(ed_f_i_combined);
      }
    }
  }
}

// original
// RBDL_DLLAPI
// void CompositeRigidBodyAlgorithm (
//   Model &model,
//   EDModel &ed_model,
//   Math::VectorNd const & q,
//   Math::MatrixNd const & q_dirs,
//   Math::MatrixNd & H,
//   std::vector<Math::MatrixNd> & H_dirs,
//   bool update_kinematics
// )
// {
//   assert (H.rows() == model.dof_count && H.cols() == model.dof_count);

//   // get number if directions
//   const size_t ndirs = q_dirs.cols();
//   ed_model.resize_directions(ndirs);

//   if (update_kinematics)
//   {
//     for (unsigned int i = 1; i < model.mBodies.size(); i++)
//     {
//       jcalc_X_lambda_S (model, i, q);
//     }
//   }

//   for (unsigned int i = 1; i < model.mBodies.size(); i++)
//   {
//     // nominal evaluation
//     model.Ic[i] = model.I[i];
//     // derivative evaluation
//     for (unsigned idir = 0; idir < ndirs; idir++)
//     {
//       ed_model.Ic[i][idir].setZero();
//     }
//   }

//   for (unsigned int i = model.mBodies.size() - 1; i > 0; i--)
//   {
//     if (model.lambda[i] != 0)
//     {
//       // nominal evaluation
//       // NOTE we require temporary spatial rbi for efficient computation
//       SpatialRigidBodyInertia temp = model.X_lambda[i].applyTranspose(model.Ic[i]);
//       model.Ic[model.lambda[i]] += temp;

//       // derivative evaluation
//       // NOTE we have to transform back S vector
//       const Math::SpatialVector imv = model.X_lambda[i].inverse().apply(model.S[i]);

//       for (unsigned idir = 0; idir < ndirs; idir++) {
//         Vector3d E_T_mr = model.X_lambda[i].E.transpose() * ed_model.Ic[i][idir].h;
//         SpatialRigidBodyInertia rbi = SpatialRigidBodyInertia (
//           0.,
//           E_T_mr,
//           model.X_lambda[i].E.transpose() *
//           Matrix3d (
//             ed_model.Ic[i][idir].Ixx, ed_model.Ic[i][idir].Iyx, ed_model.Ic[i][idir].Izx,
//             ed_model.Ic[i][idir].Iyx, ed_model.Ic[i][idir].Iyy, ed_model.Ic[i][idir].Izy,
//             ed_model.Ic[i][idir].Izx, ed_model.Ic[i][idir].Izy, ed_model.Ic[i][idir].Izz
//             ) * model.X_lambda[i].E
//           - VectorCrossMatrix(model.X_lambda[i].r) * VectorCrossMatrix (model.X_lambda[i].E.transpose() * ed_model.Ic[i][idir].h)
//           - VectorCrossMatrix (E_T_mr) * VectorCrossMatrix (model.X_lambda[i].r)
//         );

//         const Math::Vector3d w = imv.head(3)*q_dirs(model.mJoints[i].q_index, idir);
//         const Math::Vector3d v0 = imv.tail(3)*q_dirs(model.mJoints[i].q_index, idir);
//         ed_model.Ic[model.lambda[i]][idir]
//           += SpatialRigidBodyInertia (
//               0,
//               w.cross(temp.h) + temp.m * v0,
//               Matrix3d(
//                 -temp.Iyx*w[2] + temp.Izx*w[1] - temp.Iyx*w[2] + temp.Izx*w[1] + 2.*(temp.h[1]*v0[1] + temp.h[2]*v0[2]),
//                  temp.Ixx*w[2] - temp.Izx*w[0] - temp.Iyy*w[2] + temp.Izy*w[1] -     temp.h[0]*v0[1] - temp.h[1]*v0[0] ,
//                 -temp.Ixx*w[1] + temp.Iyx*w[0] - temp.Izy*w[2] + temp.Izz*w[1] -     temp.h[0]*v0[2] - temp.h[2]*v0[0] ,

//                  temp.Ixx*w[2] - temp.Iyy*w[2] + temp.Izy*w[1] - temp.Izx*w[0] -     temp.h[0]*v0[1] - temp.h[1]*v0[0] ,
//                  temp.Iyx*w[2] + temp.Iyx*w[2] - temp.Izy*w[0] - temp.Izy*w[0] + 2.*(temp.h[0]*v0[0] + temp.h[2]*v0[2]),
//                  temp.Izx*w[2] - temp.Iyx*w[1] + temp.Iyy*w[0] - temp.Izz*w[0] -     temp.h[1]*v0[2] - temp.h[2]*v0[1] ,

//                 -temp.Ixx*w[1] + temp.Iyx*w[0] - temp.Izy*w[2] + temp.Izz*w[1] -     temp.h[0]*v0[2] - temp.h[2]*v0[0] ,
//                 -temp.Iyx*w[1] + temp.Iyy*w[0] + temp.Izx*w[2] - temp.Izz*w[0] -     temp.h[1]*v0[2] - temp.h[2]*v0[1] ,
//                 -temp.Izx*w[1] + temp.Izy*w[0] - temp.Izx*w[1] + temp.Izy*w[0] + 2.*(temp.h[0]*v0[0] + temp.h[1]*v0[1])
//               )
//           )
//           + rbi
//         ;
//       }
//     }

//     unsigned int dof_index_i = model.mJoints[i].q_index;

//     if (model.mJoints[i].mDoFCount == 1 && model.mJoints[i].mJointType != JointTypeCustom)
//     {
//       // nominal evaluation
//       SpatialVector F             = model.Ic[i] * model.S[i];
//       H(dof_index_i, dof_index_i) = model.S[i].dot(F);

//       // derivative evaluation
//       for (unsigned idir = 0; idir < ndirs; idir++)
//       {
//         ed_model.F[i].col(idir) = ed_model.Ic[i][idir] * model.S[i];
//         H_dirs[idir](dof_index_i, dof_index_i) = model.S[i].dot(ed_model.F[i].col(idir));
//       }

//       unsigned int j = i;
//       unsigned int dof_index_j = dof_index_i;

//       while (model.lambda[j] != 0) {
//         // derivative evaluation
//         // TODO do not evaluate zero blocks
//         ed_model.F[i].leftCols(ndirs) = model.X_lambda[j].toMatrixTranspose()
//             * (ed_model.F[i].leftCols(ndirs) + crossf_rhs_T(F)*model.S[j]*q_dirs.row(model.mJoints[j].q_index).leftCols(ndirs));
//         // nominal evaluation
//         F = model.X_lambda[j].applyTranspose(F);

//         j = model.lambda[j];
//         dof_index_j = model.mJoints[j].q_index;

//         // nominal evaluation
//         H(dof_index_i,dof_index_j) = F.dot(model.S[j]);
//         H(dof_index_j,dof_index_i) = H(dof_index_i,dof_index_j);
//         // derivative evaluation
//         for (unsigned idir = 0; idir < ndirs; idir++) {
//           H_dirs[idir](dof_index_i, dof_index_j)
//             = ed_model.F[i].col(idir).dot(model.S[j]);
//           H_dirs[idir](dof_index_j, dof_index_i)
//             = H_dirs[idir](dof_index_i, dof_index_j);
//         }
//       }
//     } else if (model.mJoints[i].mDoFCount == 3 && model.mJoints[i].mJointType != JointTypeCustom) {
//       cerr << __FILE__ << " " << __LINE__ << ":"
//            << "Multi-DoF joint not supported." << endl;
//       abort();
//     } else if (model.mJoints[i].mJointType == JointTypeCustom) {
//       cerr << __FILE__ << " " << __LINE__ << ":"
//            << " Custom joints not supported." << endl;
//       abort();
//     } else {
//       cerr << __FILE__ << " " << __LINE__
//            << ": Unknown unsupported joint." << endl;
//       abort();
//     }
//   }
// }

// // optimized
// RBDL_DLLAPI
// void CompositeRigidBodyAlgorithm (
//   Model &model,
//   EDModel &ed_model,
//   Math::VectorNd const & q,
//   Math::MatrixNd const & q_dirs,
//   Math::MatrixNd & H,
//   std::vector<Math::MatrixNd> & H_dirs,
//   bool update_kinematics
// )
// {
//   assert (H.rows() == model.dof_count && H.cols() == model.dof_count);

//   const size_t ndirs = q_dirs.cols();
//   ed_model.resize_directions(ndirs);

//   const unsigned int num_bodies = model.mBodies.size();

//   if (update_kinematics) {
//     for (unsigned int i = 1; i < num_bodies; i++) {
//       jcalc_X_lambda_S (model, i, q);
//     }
//   }

//   // Pre-initialize nominal and derivative inertia
//   for (unsigned int i = 1; i < num_bodies; i++) {
//     model.Ic[i] = model.I[i];
//     for (unsigned idir = 0; idir < ndirs; idir++) {
//       ed_model.Ic[i][idir].setZero();
//     }
//   }

//   for (unsigned int i = num_bodies - 1; i > 0; i--) {
//     const unsigned int lambda_i = model.lambda[i];
//     const unsigned int dof_index_i = model.mJoints[i].q_index;

//     // Check unsupported joints early
//     if (model.mJoints[i].mDoFCount != 1 || model.mJoints[i].mJointType == JointTypeCustom) {
//        std::cerr << __FILE__ << ":" << __LINE__ << ": Only 1-DoF standard joints supported." << std::endl;
//        abort();
//     }

//     if (lambda_i != 0) {
//       // --- Backward Pass: Propagate Composite Inertias ---
//       SpatialTransform X_lam_i = model.X_lambda[i];
//       const Matrix3d& E_T = X_lam_i.E.transpose();
      
//       // Nominal evaluation
//       SpatialRigidBodyInertia temp = X_lam_i.applyTranspose(model.Ic[i]);
//       model.Ic[lambda_i] += temp;

//       // Pre-compute common terms for derivative evaluation
//       const Math::SpatialVector imv = X_lam_i.inverse().apply(model.S[i]);
//       const Math::Vector3d imv_w = imv.head<3>();
//       const Math::Vector3d imv_v = imv.tail<3>();

//       // Extract temp inertia properties once
//       const double tx = temp.h[0], ty = temp.h[1], tz = temp.h[2];
//       const double txx = temp.Ixx, tyy = temp.Iyy, tzz = temp.Izz;
//       const double txy = temp.Iyx, tyz = temp.Izy, txz = temp.Izx;
//       const double tm = temp.m;

//       for (unsigned idir = 0; idir < ndirs; idir++) {
//         // --- 1. Compute rbi (transformed derivative inertia) ---
//         const SpatialRigidBodyInertia& ed_I = ed_model.Ic[i][idir];
//         Vector3d E_T_mr = E_T * ed_I.h;
        
//         Matrix3d I_rot = E_T * Matrix3d(
//             ed_I.Ixx, ed_I.Iyx, ed_I.Izx,
//             ed_I.Iyx, ed_I.Iyy, ed_I.Izy,
//             ed_I.Izx, ed_I.Izy, ed_I.Izz
//         ) * X_lam_i.E;

//         SpatialRigidBodyInertia rbi(
//             0., 
//             E_T_mr, 
//             I_rot - VectorCrossMatrix(X_lam_i.r) * VectorCrossMatrix(E_T_mr) 
//                   - VectorCrossMatrix(E_T_mr) * VectorCrossMatrix(X_lam_i.r)
//         );

//         // --- 2. Compute velocity-induced inertia change ---
//         const double qd = q_dirs(dof_index_i, idir);
//         const Math::Vector3d w = imv_w * qd;
//         const Math::Vector3d v0 = imv_v * qd;
        
//         const double wx = w[0], wy = w[1], wz = w[2];
//         const double vx = v0[0], vy = v0[1], vz = v0[2];

//         // Diagonal terms optimized
//         const double diag_x = -txy*wz + txz*wy - txy*wz + txz*wy + 2.*(ty*vy + tz*vz);
//         const double diag_y =  txy*wz + txy*wz - tyz*wx - tyz*wx + 2.*(tx*vx + tz*vz);
//         const double diag_z = -txz*wy + tyz*wx - txz*wy + tyz*wx + 2.*(tx*vx + ty*vy);

//         // Off-diagonal terms optimized
//         const double xy = txx*wz - txz*wx - tyy*wz + tyz*wy - tx*vy - ty*vx;
//         const double xz = -txx*wy + txy*wx - tyz*wz + tzz*wy - tx*vz - tz*vx;
//         const double yz = -txy*wy + tyy*wx + txz*wz - tzz*wx - ty*vz - tz*vy;

//         ed_model.Ic[lambda_i][idir] += SpatialRigidBodyInertia(
//             0.,
//             w.cross(temp.h) + tm * v0,
//             Matrix3d(
//                 diag_x, xy,     xz,
//                 xy,     diag_y, yz,
//                 xz,     yz,     diag_z
//             )
//         ) + rbi;
//       }
//     }

//     // --- Compute F and H Matrix Entries ---
//     const SpatialVector& S_i = model.S[i];
//     SpatialVector F = model.Ic[i] * S_i;
//     H(dof_index_i, dof_index_i) = S_i.dot(F);

//     for (unsigned idir = 0; idir < ndirs; idir++) {
//       ed_model.F[i].col(idir) = ed_model.Ic[i][idir] * S_i;
//       H_dirs[idir](dof_index_i, dof_index_i) = S_i.dot(ed_model.F[i].col(idir));
//     }

//     unsigned int j = i;
//     unsigned int dof_index_j = dof_index_i;

//     while (model.lambda[j] != 0) {
//       const unsigned int lambda_j = model.lambda[j];
//       SpatialTransform X_lam_j = model.X_lambda[j];
//       const SpatialVector& S_j = model.S[j];
//       const unsigned int q_idx_j = model.mJoints[j].q_index;

//       // Pre-compute crossf matrix
//       Math::MatrixNd crossF_rhs_T = crossf_rhs_T(F); 
//       Math::MatrixNd X_lam_j_matT = X_lam_j.toMatrixTranspose();

//       // Derivative evaluation block
//       ed_model.F[i].leftCols(ndirs) = X_lam_j_matT * (ed_model.F[i].leftCols(ndirs) + crossF_rhs_T * S_j * q_dirs.row(q_idx_j).leftCols(ndirs));

//       // Nominal evaluation block
//       F = X_lam_j.applyTranspose(F);

//       j = lambda_j;
//       dof_index_j = model.mJoints[j].q_index;

//       // Populate H matrices
//       double h_val = F.dot(model.S[j]);
//       H(dof_index_i, dof_index_j) = h_val;
//       H(dof_index_j, dof_index_i) = h_val;

//       for (unsigned idir = 0; idir < ndirs; idir++) {
//         double hd_val = ed_model.F[i].col(idir).dot(model.S[j]);
//         H_dirs[idir](dof_index_i, dof_index_j) = hd_val;
//         H_dirs[idir](dof_index_j, dof_index_i) = hd_val;
//       }
//     }
//   }
// }

// only derivatives
RBDL_DLLAPI
void CompositeRigidBodyAlgorithm (
  Model &model,
  EDModel &ed_model,
  Math::VectorNd const & q,
  Math::MatrixNd const & q_dirs,
  std::vector<Math::MatrixNd> & H_dirs,
  bool update_kinematics
)
{
  const size_t ndirs = q_dirs.cols();
  ed_model.resize_directions(ndirs);

  const unsigned int num_bodies = model.mBodies.size();

  if (update_kinematics) {
    for (unsigned int i = 1; i < num_bodies; i++) {
      jcalc_X_lambda_S (model, i, q);
    }
  }

  // Pre-initialize nominal and derivative inertia
  for (unsigned int i = 1; i < num_bodies; i++) {
    model.Ic[i] = model.I[i];
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ed_model.Ic[i][idir].setZero();
    }
  }

  for (unsigned int i = num_bodies - 1; i > 0; i--) {
    const unsigned int lambda_i = model.lambda[i];
    const unsigned int dof_index_i = model.mJoints[i].q_index;

    // Check unsupported joints
    if (model.mJoints[i].mDoFCount != 1 || model.mJoints[i].mJointType == JointTypeCustom) {
       std::cerr << __FILE__ << ":" << __LINE__ << ": Only 1-DoF standard joints supported." << std::endl;
       abort();
    }

    if (lambda_i != 0) {
      // --- Backward Pass: Propagate Composite Inertias ---
      SpatialTransform X_lam_i = model.X_lambda[i];
      const Matrix3d& E_T = X_lam_i.E.transpose();
      
      // Nominal evaluation (Required for derivatives)
      SpatialRigidBodyInertia temp = X_lam_i.applyTranspose(model.Ic[i]);
      model.Ic[lambda_i] += temp;

      // Pre-compute common terms for derivative evaluation
      const Math::SpatialVector imv = X_lam_i.inverse().apply(model.S[i]);
      const Math::Vector3d imv_w = imv.head<3>();
      const Math::Vector3d imv_v = imv.tail<3>();

      const double tm = temp.m;
      const double tx = temp.h[0], ty = temp.h[1], tz = temp.h[2];
      const double txx = temp.Ixx, tyy = temp.Iyy, tzz = temp.Izz;
      const double txy = temp.Iyx, tyz = temp.Izy, txz = temp.Izx;

      for (unsigned idir = 0; idir < ndirs; idir++) {
        const SpatialRigidBodyInertia& ed_I = ed_model.Ic[i][idir];
        Vector3d E_T_mr = E_T * ed_I.h;
        
        Matrix3d I_rot = E_T * Matrix3d(
            ed_I.Ixx, ed_I.Iyx, ed_I.Izx,
            ed_I.Iyx, ed_I.Iyy, ed_I.Izy,
            ed_I.Izx, ed_I.Izy, ed_I.Izz
        ) * X_lam_i.E;

        SpatialRigidBodyInertia rbi(
            0., 
            E_T_mr, 
            I_rot - VectorCrossMatrix(X_lam_i.r) * VectorCrossMatrix(E_T_mr) 
                  - VectorCrossMatrix(E_T_mr) * VectorCrossMatrix(X_lam_i.r)
        );

        const double qd = q_dirs(dof_index_i, idir);
        const Math::Vector3d w = imv_w * qd;
        const Math::Vector3d v0 = imv_v * qd;
        
        const double wx = w[0], wy = w[1], wz = w[2];
        const double vx = v0[0], vy = v0[1], vz = v0[2];

        const double diag_x = -txy*wz + txz*wy - txy*wz + txz*wy + 2.*(ty*vy + tz*vz);
        const double diag_y =  txy*wz + txy*wz - tyz*wx - tyz*wx + 2.*(tx*vx + tz*vz);
        const double diag_z = -txz*wy + tyz*wx - txz*wy + tyz*wx + 2.*(tx*vx + ty*vy);

        const double xy = txx*wz - txz*wx - tyy*wz + tyz*wy - tx*vy - ty*vx;
        const double xz = -txx*wy + txy*wx - tyz*wz + tzz*wy - tx*vz - tz*vx;
        const double yz = -txy*wy + tyy*wx + txz*wz - tzz*wx - ty*vz - tz*vy;

        ed_model.Ic[lambda_i][idir] += SpatialRigidBodyInertia(
            0.,
            w.cross(temp.h) + tm * v0,
            Matrix3d(
                diag_x, xy,     xz,
                xy,     diag_y, yz,
                xz,     yz,     diag_z
            )
        ) + rbi;
      }
    }

    // --- Compute F and H_dirs Matrix Entries ---
    const SpatialVector& S_i = model.S[i];
    SpatialVector F = model.Ic[i] * S_i; 

    for (unsigned idir = 0; idir < ndirs; idir++) {
      ed_model.F[i].col(idir) = ed_model.Ic[i][idir] * S_i;
      H_dirs[idir](dof_index_i, dof_index_i) = S_i.dot(ed_model.F[i].col(idir));
    }

    unsigned int j = i;
    while (model.lambda[j] != 0) {
      const unsigned int lambda_j = model.lambda[j];
      SpatialTransform X_lam_j = model.X_lambda[j];
      const SpatialVector& S_j = model.S[j];
      const unsigned int q_idx_j = model.mJoints[j].q_index;

      // --- OPTIMIZED DERIVATIVE PASS ---
      // 1. Efficient Spatial Force Cross Product: F_cross = S_j x* F
      SpatialVector F_cross = crossf(S_j, F);

      for (unsigned idir = 0; idir < ndirs; idir++) {
        // 2. Linear combination: F_dir = F_dir + (S_j x* F) * q_dir
        ed_model.F[i].col(idir) += F_cross * q_dirs(q_idx_j, idir);
        
        // 3. Coordinate transform: F_dir = X^T * F_dir
        // applyTranspose is O(n) vs O(n^2) for toMatrixTranspose() * vector
        ed_model.F[i].col(idir) = X_lam_j.applyTranspose(ed_model.F[i].col(idir));
      }

      // 4. Update nominal force (Required for next cross product)
      F = X_lam_j.applyTranspose(F);

      j = lambda_j;
      const unsigned int dof_index_j = model.mJoints[j].q_index;

      // 5. Populate H_dirs entries (H removed)
      for (unsigned idir = 0; idir < ndirs; idir++) {
        double hd_val = ed_model.F[i].col(idir).dot(model.S[j]);
        H_dirs[idir](dof_index_i, dof_index_j) = hd_val;
        H_dirs[idir](dof_index_j, dof_index_i) = hd_val;
      }
    }
  }
}

// // optimized derivatives
// RBDL_DLLAPI
// void CompositeRigidBodyAlgorithm (
//   Model &model,
//   EDModel &ed_model,
//   Math::VectorNd const & q,
//   Math::MatrixNd const & q_dirs,
//   std::vector<Math::MatrixNd> & H_dirs,
//   bool update_kinematics
// )
// {
//   const size_t ndirs = q_dirs.cols();
//   ed_model.resize_directions(ndirs);

//   const unsigned int num_bodies = model.mBodies.size();

//   if (update_kinematics) {
//     for (unsigned int i = 1; i < num_bodies; i++) {
//       jcalc_X_lambda_S (model, i, q);
//     }
//   }

//   // Pre-initialize nominal and derivative inertia
//   for (unsigned int i = 1; i < num_bodies; i++) {
//     model.Ic[i] = model.I[i];
//     for (unsigned idir = 0; idir < ndirs; idir++) {
//       ed_model.Ic[i][idir].setZero();
//     }
//   }

//   for (unsigned int i = num_bodies - 1; i > 0; i--) {
//     const unsigned int lambda_i = model.lambda[i];
//     const unsigned int dof_index_i = model.mJoints[i].q_index;

//     // Check unsupported joints
//     if (model.mJoints[i].mDoFCount != 1 || model.mJoints[i].mJointType == JointTypeCustom) {
//        std::cerr << __FILE__ << ":" << __LINE__ << ": Only 1-DoF standard joints supported." << std::endl;
//        abort();
//     }

//     if (lambda_i != 0) {
//       // --- Backward Pass: Propagate Composite Inertias ---
//       SpatialTransform X_lam_i = model.X_lambda[i];
//       const Matrix3d E_T = X_lam_i.E.transpose();
//       const Vector3d& r_i = X_lam_i.r;
      
//       // Nominal evaluation (Required for derivatives)
//       SpatialRigidBodyInertia temp = X_lam_i.applyTranspose(model.Ic[i]);
//       model.Ic[lambda_i] += temp;

//       // Pre-compute common terms for derivative evaluation
//       const Math::SpatialVector imv = X_lam_i.inverse().apply(model.S[i]);
//       const Math::Vector3d imv_w = imv.head<3>();
//       const Math::Vector3d imv_v = imv.tail<3>();

//       const double tm = temp.m;
//       const double tx = temp.h[0], ty = temp.h[1], tz = temp.h[2];
//       const double txx = temp.Ixx, tyy = temp.Iyy, tzz = temp.Izz;
//       const double txy = temp.Iyx, tyz = temp.Izy, txz = temp.Izx;

//       for (unsigned idir = 0; idir < ndirs; idir++) {
//         const SpatialRigidBodyInertia& ed_I = ed_model.Ic[i][idir];
//         Vector3d E_T_mr = E_T * ed_I.h;
        
//         Matrix3d I_mat;
//         I_mat << ed_I.Ixx, ed_I.Iyx, ed_I.Izx,
//                  ed_I.Iyx, ed_I.Iyy, ed_I.Izy,
//                  ed_I.Izx, ed_I.Izy, ed_I.Izz;

//         Matrix3d I_rot = E_T * I_mat * X_lam_i.E;

//         // Analytic simplification for symmetric cross-product matrices
//         double r_dot_mr = r_i.dot(E_T_mr);
//         Matrix3d sym_cross = r_i * E_T_mr.transpose() + E_T_mr * r_i.transpose();
//         sym_cross(0,0) -= 2.0 * r_dot_mr;
//         sym_cross(1,1) -= 2.0 * r_dot_mr;
//         sym_cross(2,2) -= 2.0 * r_dot_mr;

//         SpatialRigidBodyInertia rbi(0., E_T_mr, I_rot - sym_cross);

//         const double qd = q_dirs(dof_index_i, idir);
//         const Math::Vector3d w = imv_w * qd;
//         const Math::Vector3d v0 = imv_v * qd;
        
//         const double wx = w[0], wy = w[1], wz = w[2];
//         const double vx = v0[0], vy = v0[1], vz = v0[2];

//         // Algebraic reductions
//         const double diag_x = 2.0 * (-txy*wz + txz*wy + ty*vy + tz*vz);
//         const double diag_y = 2.0 * ( txy*wz - tyz*wx + tx*vx + tz*vz);
//         const double diag_z = 2.0 * (-txz*wy + tyz*wx + tx*vx + ty*vy);

//         const double xy = txx*wz - txz*wx - tyy*wz + tyz*wy - tx*vy - ty*vx;
//         const double xz = -txx*wy + txy*wx - tyz*wz + tzz*wy - tx*vz - tz*vx;
//         const double yz = -txy*wy + tyy*wx + txz*wz - tzz*wx - ty*vz - tz*vy;

//         ed_model.Ic[lambda_i][idir] += SpatialRigidBodyInertia(
//             0.,
//             w.cross(temp.h) + tm * v0,
//             Matrix3d(
//                 diag_x, xy,     xz,
//                 xy,     diag_y, yz,
//                 xz,     yz,     diag_z
//             )
//         ) + rbi;
//       }
//     }

//     // --- Compute F and H_dirs Matrix Entries ---
//     const SpatialVector& S_i = model.S[i];
//     SpatialVector F = model.Ic[i] * S_i; 

//     for (unsigned idir = 0; idir < ndirs; idir++) {
//       ed_model.F[i].col(idir) = ed_model.Ic[i][idir] * S_i;
//       H_dirs[idir](dof_index_i, dof_index_i) = S_i.dot(ed_model.F[i].col(idir));
//     }

//     unsigned int j = i;
//     while (model.lambda[j] != 0) {
//       const unsigned int lambda_j = model.lambda[j];
//       SpatialTransform X_lam_j = model.X_lambda[j];
      
//       // Keep S_j local for the F_cross calculation
//       const SpatialVector& S_j_old = model.S[j]; 
//       const unsigned int q_idx_j = model.mJoints[j].q_index;

//       SpatialVector F_cross = crossf(S_j_old, F);

//       // Vectorized Rank-1 Update
//       ed_model.F[i].noalias() += F_cross * q_dirs.row(q_idx_j);

//       for (unsigned idir = 0; idir < ndirs; idir++) {
//         ed_model.F[i].col(idir) = X_lam_j.applyTranspose(ed_model.F[i].col(idir));
//       }

//       F = X_lam_j.applyTranspose(F);

//       // Step up the tree
//       j = lambda_j;
//       const unsigned int dof_index_j = model.mJoints[j].q_index;

//       // FIX: Use model.S[j] here so we're evaluating with the newly updated 'j'
//       Eigen::RowVectorXd hd_vals = model.S[j].transpose() * ed_model.F[i];
//       for (unsigned idir = 0; idir < ndirs; idir++) {
//         H_dirs[idir](dof_index_i, dof_index_j) = hd_vals(idir);
//         H_dirs[idir](dof_index_j, dof_index_i) = hd_vals(idir);
//       }
//     }
//   }
// }

// -----------------------------------------------------------------------------
} // ED
// -----------------------------------------------------------------------------
} // RigidBodyDynamics
// -----------------------------------------------------------------------------
