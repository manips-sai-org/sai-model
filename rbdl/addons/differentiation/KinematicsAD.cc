/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <limits>
#include <cstring>
#include <assert.h>

#include "rbdl/rbdl_mathutils.h"
#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"

#include "JointAD.h"
#include "ModelAD.h"
#include "KinematicsAD.h"
#include "SpatialAlgebraOperatorsAD.h"

#include "rbdl_mathutilsAD.h"
#include "rbdl_mathutilsFD.h"

using std::cerr;
using std::endl;
using std::vector;

using namespace RigidBodyDynamics::Math;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI Vector3d CalcBodyToBaseCoordinates (
    Model &model,
    ADModel &ad_model,
    const VectorNd &Q,
    const MatrixNd &Q_dirs,
    unsigned int body_id,
    const Vector3d &point_body_coordinates,
    Math::MatrixNd & ad_body_to_base_coordinates,
    bool update_kinematics
    ) {
  unsigned ndirs = Q_dirs.cols();

  vector<Matrix3d> ad_body_rotation (ndirs, Matrix3d::Zero());
  vector<Vector3d> ad_body_position (ndirs, Vector3d::Zero());

  if (update_kinematics) {
    // nominal + derivative evaluation
    AD::UpdateKinematicsCustom (
          model, ad_model,
          &Q, &Q_dirs,
          NULL, NULL,
          NULL, NULL
          );
  }

  if (body_id >= model.fixed_body_discriminator) {
    cerr << __FILE__ << " " << __LINE__
         << ": Fixed bodies not yet supported!" << endl;
    abort();
  }

  // derivative evaluation
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_body_rotation[idir] = ad_model.X_base[body_id][idir].E.transpose();
    ad_body_position[idir] = ad_model.X_base[body_id][idir].r;

    // NOTE point_body_coordinates is constant, no derivative needed!
    ad_body_to_base_coordinates.col(idir) = ad_body_position[idir]
        + ad_body_rotation[idir] * point_body_coordinates;
  }

  // nominal evaluation
  Matrix3d body_rotation = model.X_base[body_id].E.transpose();
  Vector3d body_position = model.X_base[body_id].r;
  return body_position + body_rotation * point_body_coordinates;
}


RBDL_DLLAPI Vector3d CalcBaseToBodyCoordinates (
    Model & model,
    ADModel & ad_model,
    VectorNd const & q,
    MatrixNd const & q_dirs,
    unsigned body_id,
    Vector3d const & base_pt_pos,
    MatrixNd const & base_pt_pos_dirs,
    MatrixNd & ad_base_to_body_coordinates,
    bool update_kinematics) {
  unsigned ndirs = q_dirs.cols();
  assert(ndirs == base_pt_pos_dirs.cols());
  assert(ndirs == ad_base_to_body_coordinates.cols());
  assert(3     == base_pt_pos_dirs.rows());
  assert(3     == ad_base_to_body_coordinates.rows());

  if (update_kinematics) {
    UpdateKinematicsCustom (model, ad_model, &q, &q_dirs, 0, 0, 0, 0);
  }

  if (body_id >= model.fixed_body_discriminator) {
    cerr << __FILE__ << " " << __LINE__
         << ": Fixed bodies not yet supported!" << endl;
    abort();
  }

  vector<Matrix3d> ad_body_rotation (ndirs, Matrix3d::Zero());
  MatrixNd ad_body_position (3, ndirs);

  // derivative code
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_body_rotation[idir] = ad_model.X_base[body_id][idir].E;
    ad_body_position.col(idir) = ad_model.X_base[body_id][idir].r;
  }
  // nominal code
  Matrix3d body_rotation = model.X_base[body_id].E;
  Vector3d body_position = model.X_base[body_id].r;

  // derivative code
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_base_to_body_coordinates.col(idir) =
        ad_body_rotation[idir] * (base_pt_pos - body_position)
        + body_rotation * (base_pt_pos_dirs.col(idir)
                           - ad_body_position.col(idir));
  }
  // nominal code
  return body_rotation * (base_pt_pos - body_position);
}


RBDL_DLLAPI void UpdateKinematicsCustom (
    Model &model,
    ADModel & ad_model,
    const VectorNd * q,
    const MatrixNd * q_dirs,
    const VectorNd * qdot,
    const MatrixNd * qdot_dirs,
    const VectorNd * qddot,
    const MatrixNd * qddot_dirs) {
  LOG << "-------- " << __func__ << " --------" << std::endl;
  unsigned int i;
  unsigned int ndirs = 0;

  if (q) {
    ndirs = q_dirs->cols();
    assert(!qdot || qdot_dirs->cols() == ndirs);
    assert(!qddot || qddot_dirs->cols() == ndirs);
  } else if(qdot) {
    ndirs = qdot_dirs->cols();
    assert(!qddot || qddot_dirs->cols() == ndirs);
  } else if (qddot) {
    ndirs = qddot_dirs->cols();
  }

  ad_model.resize_directions(ndirs);

  if (q) {
    for (i = 1; i < model.mBodies.size(); i++) {
      unsigned int lambda = model.lambda[i];

      VectorNd QDot_zero (model.q_size);
      QDot_zero.setZero();
      MatrixNd QDot_zero_dirs (model.q_size, ndirs);
      QDot_zero_dirs.setZero();

      // Derivative evaluation and nominal evaluation
      jcalc (model, ad_model, i, *q, *q_dirs, QDot_zero, QDot_zero_dirs);

//      // derivative evaluation
//      for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
//        // NOTE: X_T is a constant model dependent transformation
//        ad_model.X_lambda[i][idirs] = ad_model.X_J[i][idirs] * model.X_T[i].toMatrix();
//      }
//      // nominal evaluation
//      model.X_lambda[i] = model.X_J[i] * model.X_T[i];
      mulSTST(ndirs,
              model.X_J[i], ad_model.X_J[i],
              model.X_T[i],
              model.X_lambda[i], ad_model.X_lambda[i]);

      if (lambda != 0) {
//        // derivative evaluation
//        for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
//          ad_model.X_base[i][idirs] =
//              ad_model.X_lambda[i][idirs] * model.X_base[lambda].toMatrix()
//              + model.X_lambda[i].toMatrix() * ad_model.X_base[lambda][idirs];
//        }
//        // nominal evaluation
//        model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];

        mulSTST(ndirs,
                model.X_lambda[i], ad_model.X_lambda[i],
                model.X_base[lambda], ad_model.X_base[lambda],
                model.X_base[i], ad_model.X_base[i]);

      } else {
        // derivative evaluation
        for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
          ad_model.X_base[i][idirs] = ad_model.X_lambda[i][idirs];
        }
        // nominal evaluation
        model.X_base[i] = model.X_lambda[i];
      }
    }
  }

  if (qdot) {
    for (i = 1; i < model.mBodies.size(); i++) {
      unsigned int lambda = model.lambda[i];

      // Derivative evaluation and nominal evaluation
      jcalc (model, ad_model, i, *q, *q_dirs, *qdot, *qdot_dirs);

      if (lambda != 0) {
        applySTSV(ndirs,
                  model.X_lambda[i], ad_model.X_lambda[i],
                  model.v[lambda], ad_model.v[lambda],
                  model.v[i], ad_model.v[i]);

        for (unsigned idir = 0; idir < ndirs; idir++) {
          ad_model.v[i][idir] += ad_model.v_J[i][idir];
        }
        model.v[i] += model.v_J[i];
      } else {
        // derivative evaluation
        for (unsigned idir = 0; idir < ndirs; idir++) {
          ad_model.v[i][idir] = ad_model.v_J[i][idir];
        }
        // nominal evaluation
        model.v[i] = model.v_J[i];
      }
      // derivative evaluation
      for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
        ad_model.c[i][idirs] = ad_model.c_J[i][idirs]
            + AD::crossm(
              model.v[i], ad_model.v[i][idirs],
              model.v_J[i], ad_model.v_J[i][idirs]);
      }
      // nominal evaluation
      model.c[i] = model.c_J[i] + Math::crossm(model.v[i],model.v_J[i]);
    }
  }

  if (qddot) {
    for (i = 1; i < model.mBodies.size(); i++) {
      unsigned int q_index = model.mJoints[i].q_index;
      unsigned int lambda = model.lambda[i];

      if (lambda != 0) {
        applySTSV(ndirs,
                  model.X_lambda[i], ad_model.X_lambda[i],
                  model.a[lambda], ad_model.a[lambda],
                  model.a[i], ad_model.a[i]);
        for (unsigned idir = 0; idir < ndirs; idir++) {
          ad_model.a[i][idir] += ad_model.c[i][idir];
        }
        model.a[i] += model.c[i];
      } else {
        // derivative evaluation
        for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
          ad_model.a[i][idirs] = ad_model.c[i][idirs];
        }
        // nominal evaluation
        model.a[i] = model.c[i];
      }

      if (model.mJoints[i].mDoFCount == 3) {
        cerr << "Multi-DoF not supported." << endl;
        abort();
      } else {
        // derivative evaluation
        for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
          ad_model.a[i][idirs] = ad_model.a[i][idirs]
              + model.S[i] * (*qddot_dirs)(q_index, idirs);
        }
        // nominal evaluation
        model.a[i] = model.a[i] + model.S[i] * (*qddot)[q_index];
      }
    }
  }
}

RBDL_DLLAPI void UpdateKinematics (
    Model &model,
    ADModel &ad_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs,
    const VectorNd &qddot,
    const MatrixNd &qddot_dirs) {
  LOG << "-------- " << __func__ << " --------" << std::endl;

  unsigned int i;
  unsigned int ndirs = q_dirs.cols();

  ad_model.resize_directions(ndirs);

  // derivative evaluation
  for (unsigned int idir = 0; idir < ndirs; ++idir) {
    ad_model.a[0][idir].setZero();
  }
  // nominal evaluation
  model.a[0].setZero();
  //model.a[0] = spatial_gravity;

  for (i = 1; i < model.mBodies.size(); i++) {
    unsigned int q_index = model.mJoints[i].q_index;

    Joint joint = model.mJoints[i];
    unsigned int lambda = model.lambda[i];

    // nominal + derivative evaluation
    jcalc (model, ad_model, i, q, q_dirs, qdot, qdot_dirs);

    mulSTST(ndirs,
            model.X_J[i], ad_model.X_J[i],
            model.X_T[i],
            model.X_lambda[i], ad_model.X_lambda[i]);

    if (lambda != 0) {
      mulSTST(ndirs,
              model.X_lambda[i], ad_model.X_lambda[i],
              model.X_base[lambda], ad_model.X_base[lambda],
              model.X_base[i], ad_model.X_base[i]);

      applySTSV(ndirs,
                model.X_lambda[i], ad_model.X_lambda[i],
                model.v[lambda], ad_model.v[lambda],
                model.v[i], ad_model.v[i]);
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.v[i][idir] += ad_model.v_J[i][idir];
      }
      model.v[i] += model.v_J[i];

    } else {
      // derivative evaluation
      for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
        ad_model.X_base[i][idirs] = ad_model.X_lambda[i][idirs];
        ad_model.v[i][idirs] = ad_model.v_J[i][idirs];
      }
      // nominal evaluation
      model.X_base[i] = model.X_lambda[i];
      model.v[i] = model.v_J[i];
    }

    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; ++idir) {
      ad_model.c[i][idir] = ad_model.c_J[i][idir]
          + AD::crossm(
            model.v[i], ad_model.v[i][idir],
            model.v_J[i], ad_model.v_J[i][idir]
            );
    }
    // nominal evaluation
    model.c[i] = model.c_J[i] + Math::crossm(model.v[i],model.v_J[i]);

    applySTSV(ndirs,
              model.X_lambda[i], ad_model.X_lambda[i],
              model.a[lambda], ad_model.a[lambda],
              model.a[i], ad_model.a[i]);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_model.a[i][idir] += ad_model.c[i][idir];
    }
    model.a[i] += model.c[i];

    if (model.mJoints[i].mDoFCount == 3
        || model.mJoints[i].mJointType == JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__
           << ": Multi-DoF and custom joints not supported!" << endl;
      abort();
    } else {
      // derivative evaluation
      for (unsigned int idir = 0; idir < ndirs; ++idir) {
        ad_model.a[i][idir] = ad_model.a[i][idir]
            + model.S[i] * qddot_dirs(q_index, idir);
      }
      // nominal evaluation
      model.a[i] = model.a[i] + model.S[i] * qddot[q_index];
    }
  }

  for (i = 1; i < model.mBodies.size(); i++) {
    // derivative evaluation
    for (unsigned int idirs = 0; idirs < ndirs; ++idirs) {
      LOG << "ad_a[" << i << "][" << idirs << "] = "
          << ad_model.a[i][idirs].transpose() << std::endl;
    }
    // nominal evaluation
    LOG << "a[" << i << "] = " << model.a[i].transpose() << std::endl;
  }
}


RBDL_DLLAPI Matrix3d CalcBodyWorldOrientation (
    Model & model,
    ADModel & ad_model,
    VectorNd const & q,
    MatrixNd const & q_dirs,
    const unsigned int body_id,
    vector<Matrix3d> & ad_derivative,
    bool update_kinematics) {
  unsigned ndirs = q_dirs.cols();
  assert(ad_derivative.size() == ndirs);

  // update the Kinematics if necessary
  if (update_kinematics) {
    UpdateKinematicsCustom (model, ad_model, &q, &q_dirs, 0, 0, 0, 0);
  }

  if (body_id >= model.fixed_body_discriminator) {
    std::cerr << "Fixed bodies not yet supported!" << std::endl;
    abort();
  }

  for (unsigned int idir = 0; idir < ndirs; idir++) {
    ad_derivative[idir] = ad_model.X_base[body_id][idir].E;
  }
  // nominal evaluation
  return model.X_base[body_id].E;
}


RBDL_DLLAPI Vector3d CalcPointVelocity (
    Model & model,
    ADModel & ad_model,
    VectorNd const & q,
    MatrixNd const & q_dirs,
    VectorNd const & qdot,
    MatrixNd const & qdot_dirs,
    unsigned int body_id,
    Vector3d const & point_position,
    MatrixNd & ad_point_velocity,
    bool update_kinematics) {
  unsigned ndirs = q_dirs.cols();

  assert(qdot_dirs.cols() == ndirs);
  assert (model.IsBodyId(body_id));
  assert (model.q_size == q.rows());
  assert (model.qdot_size == qdot.rows());

  // Reset the velocity of the root body
  // derivative code
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_model.v[0][idir].setZero();
  }
  // nominal code
  model.v[0].setZero();


  // update the Kinematics with zero acceleration
  if (update_kinematics) {
    UpdateKinematicsCustom (model, ad_model, &q, &q_dirs, &qdot, &qdot_dirs,
                            0, 0);
  }

  unsigned int reference_body_id = body_id;
  Vector3d reference_point = point_position;

  MatrixNd ad_base_coords(3, ndirs);
  MatrixNd ad_reference_point(3, ndirs);

  if (model.IsFixedBodyId(body_id)) {
    unsigned fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
    // nominal and derivative code
    Vector3d base_coords = CalcBodyToBaseCoordinates (model, ad_model, q, q_dirs,
                                                      body_id, point_position, ad_base_coords, false);

    // nominal and derivative code
    reference_point = CalcBaseToBodyCoordinates (model, ad_model, q, q_dirs,
                                                 reference_body_id, base_coords, ad_base_coords, ad_reference_point,
                                                 false);
  }

  // derivative and nominal code
  vector<Matrix3d> ad_bwo(ndirs, Matrix3d::Zero());
  Matrix3d bwo = CalcBodyWorldOrientation(model, ad_model, q, q_dirs,
                                          reference_body_id, ad_bwo, false);

  // derivative code
  vector<SpatialMatrix> ad_st(ndirs);
  for (unsigned idir = 0; idir < ndirs; idir++) {
    Matrix3d ad_Erx = ad_bwo[idir].transpose() * Matrix3d(
          0., -reference_point[2],  reference_point[1],
        reference_point[2],                  0., -reference_point[0],
        -reference_point[1],  reference_point[0],                  0.);
    ad_st[idir].block<3,3>(0, 0) = ad_bwo[idir].transpose();
    ad_st[idir].block<3,3>(0, 3) = Matrix3dZero;
    ad_st[idir].block<3,3>(3, 0) = -ad_Erx;
    ad_st[idir].block<3,3>(3, 3) = ad_bwo[idir].transpose();
  }

  // nominal code
  SpatialTransform st(bwo.transpose(), reference_point);

  vector<SpatialVector> ad_point_spatial_velocity(ndirs);
  // derivative code
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_point_spatial_velocity[idir] = ad_st[idir] * model.v[reference_body_id]
        + st.apply(ad_model.v[reference_body_id][idir]);
  }
  // nominal code
  SpatialVector point_spatial_velocity =
      st.apply(model.v[reference_body_id]);

  // derivative code
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_point_velocity.col(idir) = ad_point_spatial_velocity[idir].segment(3, 3);
  }
  // nominal code
  return Vector3d (
      point_spatial_velocity[3],
      point_spatial_velocity[4],
      point_spatial_velocity[5]
      );
}

RBDL_DLLAPI SpatialVector CalcPointVelocity6D(
    Model &model,
    ADModel &ad_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs,
    unsigned int body_id,
    const Vector3d &point_position,
    vector<SpatialVector> &pv6d_dirs,
    bool update_kinematics) {
  unsigned ndirs = q_dirs.cols();
  assert(ndirs == static_cast<unsigned>(qdot_dirs.cols()));
  assert(ndirs == static_cast<unsigned>(pv6d_dirs.size()));

  assert (model.IsBodyId(body_id));
  assert (model.q_size == q.size());
  assert (model.qdot_size == qdot.size());

  // Reset the velocity of the root body
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_model.v[0][idir].setZero();
  }
  model.v[0].setZero();

  // update the Kinematics with zero acceleration
  if (update_kinematics) {
    UpdateKinematicsCustom (model, ad_model,
                            &q, &q_dirs,
                            &qdot, &qdot_dirs,
                            NULL, NULL);
  }

  unsigned int reference_body_id = body_id;
  MatrixNd reference_point_dirs = MatrixNd::Zero(3, ndirs);
  Vector3d reference_point = point_position;

  if (model.IsFixedBodyId(body_id)) {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;

    MatrixNd base_coords_dirs(3, ndirs);
    Vector3d base_coords = CalcBodyToBaseCoordinates(
          model, ad_model, q, q_dirs, body_id,
          point_position, base_coords_dirs, false);

    reference_point = CalcBaseToBodyCoordinates(
          model, ad_model, q, q_dirs, reference_body_id,
          base_coords, base_coords_dirs, reference_point_dirs, false);
  }


  vector<Matrix3d> E_dirs(ndirs);
  Matrix3d E = CalcBodyWorldOrientation (model, ad_model, q, q_dirs,
      reference_body_id, E_dirs, false);
  vector<SpatialTransform> st_dirs(ndirs);
  for (unsigned idir = 0; idir < ndirs; idir++) {
    st_dirs[idir] = SpatialTransform(E_dirs[idir].transpose(),
                                     reference_point_dirs.col(idir));
  }
  SpatialTransform st(E.transpose(), reference_point);

  SpatialVector pv6d;
  applySTSV(ndirs,
            st, st_dirs,
            model.v[reference_body_id], ad_model.v[reference_body_id],
            pv6d, pv6d_dirs);
  return pv6d;
}

RBDL_DLLAPI Vector3d CalcPointAcceleration (
    Model &model,
    ADModel &ad_model,
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    const Math::VectorNd &qdot,
    const Math::MatrixNd &qdot_dirs,
    const Math::VectorNd &qddot,
    const Math::MatrixNd &qddot_dirs,
    unsigned int body_id,
    const Math::Vector3d &point_position,
    Math::MatrixNd &ad_derivative,
    bool update_kinematics
) {
  LOG << "-------- " << __func__ << " --------" << std::endl;

  int const ndirs = q_dirs.cols();
//  assert(ndirs == ad_derivative.cols());
  assert(3     == ad_derivative.rows());

  ad_model.resize_directions(ndirs);

  // Reset the velocity of the root body
  model.v[0].setZero();
  model.a[0].setZero();

  if (update_kinematics) {
    // derivative evaluation
    UpdateKinematics (model, ad_model, q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs);
    // nominal evaluation
    // NOTE: Kinematics are already updated in ad_UpdateKinematics
    //UpdateKinematics (model, q, qdot, qddot);
  }

  LOG << std::endl;

  unsigned int reference_body_id = body_id;
  Vector3d reference_point = point_position;

  if (model.IsFixedBodyId(body_id)) {
    std::cerr << "Fixed bodies not yet supported!" << std::endl;
    abort();
  }

  vector<Matrix3d> ad_E(ndirs);
  Matrix3d E = CalcBodyWorldOrientation(
    model, ad_model, q, q_dirs, reference_body_id, ad_E, false
  );

  // derivative evaluation
  vector<SpatialMatrix> ad_p_X_i(ndirs);
  for (int idir = 0; idir < ndirs; idir++) {
    Matrix3d ad_ETrx = ad_E[idir].transpose() * Matrix3d(
          0., -reference_point[2],  reference_point[1],
        reference_point[2],                  0., -reference_point[0],
        -reference_point[1],  reference_point[0],                  0.);
    // + E.transpose() * Matrix3dZero;
    ad_p_X_i[idir].block<3,3>(0, 0) = ad_E[idir].transpose();
    ad_p_X_i[idir].block<3,3>(0, 3) = Matrix3dZero;
    ad_p_X_i[idir].block<3,3>(3, 0) = -ad_ETrx;
    ad_p_X_i[idir].block<3,3>(3, 3) = ad_E[idir].transpose();
  }
  // nominal evaluation
  SpatialTransform p_X_i (E.transpose(), reference_point);

  // derivative evaluation
  vector<SpatialVector> ad_p_v_i(ndirs);
  for (int idir = 0; idir < ndirs; idir++) {
    ad_p_v_i[idir] = ad_p_X_i[idir] * model.v[reference_body_id]
        + p_X_i.apply(ad_model.v[reference_body_id][idir]);
  }
  // nominal evaluation
  SpatialVector p_v_i = p_X_i.apply(model.v[reference_body_id]);

  Vector3d p_v_i_0t2(p_v_i[0], p_v_i[1], p_v_i[2]);
  Vector3d p_v_i_3t5(p_v_i[3], p_v_i[4], p_v_i[5]);
  // derivative evaluation
  vector<Vector3d> ad_a_dash(ndirs);
  for (int idir = 0; idir < ndirs; idir++) {
    Vector3d ad_p_v_i_0t2_idir = ad_p_v_i[idir].block<3,1>(0, 0);
    Vector3d ad_p_v_i_3t5_idir = ad_p_v_i[idir].block<3,1>(3, 0);
    ad_a_dash[idir] = p_v_i_0t2.cross(ad_p_v_i_3t5_idir)
        + ad_p_v_i_0t2_idir.cross(p_v_i_3t5);
  }
  // nominal evaluation
  Vector3d a_dash = p_v_i_0t2.cross(p_v_i_3t5);

  // derivative evaluation
  vector<SpatialVector> ad_p_a_i(ndirs);
  for (int idir = 0; idir < ndirs; idir++) {
    ad_p_a_i[idir] = ad_p_X_i[idir] * model.a[reference_body_id]
        + p_X_i.apply(ad_model.a[reference_body_id][idir]);
  }
  // nominal evaluation
  SpatialVector p_a_i = p_X_i.apply(model.a[reference_body_id]);

  // derivative evaluation
  for (int idir = 0; idir < ndirs; idir++) {
    ad_derivative(0, idir) = ad_p_a_i[idir][3] + ad_a_dash[idir][0];
    ad_derivative(1, idir) = ad_p_a_i[idir][4] + ad_a_dash[idir][1];
    ad_derivative(2, idir) = ad_p_a_i[idir][5] + ad_a_dash[idir][2];
  }
  // nominal evaluation
  return Vector3d (
      p_a_i[3] + a_dash[0],
      p_a_i[4] + a_dash[1],
      p_a_i[5] + a_dash[2]);
}

RBDL_DLLAPI
void CalcPointJacobian (
    Model &model,
    ADModel &ad_model,
    VectorNd const &Q,
    MatrixNd const &Q_dirs,
    unsigned body_id,
    Vector3d const &point_position,
    MatrixNd &G,
    vector<MatrixNd> &G_dirs,
    bool update_kinematics
    ) {
  LOG << "-------- " << __func__ << " --------" << std::endl;

  unsigned int ndirs = Q_dirs.cols();
  ad_model.resize_directions(ndirs);

  // update the Kinematics if necessary
  if (update_kinematics) {
    // derivative evaluation
    UpdateKinematicsCustom (
          model, ad_model,
          &Q, &Q_dirs,
          NULL, NULL,
          NULL, NULL
          );
    // NOTE nominal evaluatioN: nominal kinematics are already updated in the
    // AD version of UpdateKinematicsCustom
  }

  // NOTE we split the evaluation of the derivatives. Therefore, we first
  //      derive CalcBodyToBaseCoordinates then the spatial transform!
  MatrixNd ad_r(3, ndirs);
  Vector3d r = AD::CalcBodyToBaseCoordinates (
        model, ad_model, Q, Q_dirs, body_id, point_position, ad_r, false
        );

  // derivative evaluation
  // vector<SpatialTransform> ad_model.X_0 (ndirs, SpatialTransform::Zero());
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_model.X_0[idir].r = ad_r.col(idir);
  }
  // nominal evaluation
  // NOTE vector r is already evaluated in derivative evaluation
  SpatialTransform point_trans = SpatialTransform (Matrix3d::Identity(), r);

  // derivative evaluation
  for (unsigned int idirs = 0; idirs < ndirs; idirs++) {
    assert (G_dirs[idirs].rows() == 3);
    assert (G_dirs[idirs].cols() == model.qdot_size );
  }
  // nominal evaluation
  assert (G.rows() == 3 && G.cols() == model.qdot_size );

  unsigned int reference_body_id = body_id;

  if (model.IsFixedBodyId(body_id)) {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
  }

  unsigned int j = reference_body_id;

  // e[j] is set to 1 if joint j contributes to the Jacobian that we are
  // computing. For all other joints the column will be zero.
  while (j != 0) {
    unsigned int q_index = model.mJoints[j].q_index;

    if(model.mJoints[j].mJointType != JointTypeCustom) {
      if (model.mJoints[j].mDoFCount == 3) {
        cerr << __FILE__ << " " << __LINE__ << ": "
             << "Multi-DoF joints not supported." << endl;
        abort();
        G.block(0, q_index, 3, 3) =
          ((point_trans
            * model.X_base[j].inverse()).toMatrix()
           * model.multdof3_S[j]).block(3,0,3,3);
      } else {
        SpatialTransform inv_X_base_j;
        vector<SpatialTransform> ad_inv_X_base_j(ndirs);

        inverse(ndirs,
                model.X_base[j], ad_model.X_base[j],
                inv_X_base_j, ad_inv_X_base_j);

        SpatialVector v;
        vector<SpatialVector> ad_v(ndirs);

        applySTSV(ndirs,
                  inv_X_base_j, ad_inv_X_base_j,
                  model.S[j],
                  v, ad_v);

        SpatialVector ptv;
        vector<SpatialVector> ad_ptv(ndirs);
        applySTSV(ndirs,
                  point_trans, ad_model.X_0,
                  v, ad_v,
                  ptv, ad_ptv);

        for (unsigned idir = 0; idir < ndirs; idir++) {
          G_dirs[idir].block<3,1>(0, q_index) = ad_ptv[idir].segment<3>(3);
        }
        G.block<3,1>(0, q_index) = point_trans.apply(v).segment<3>(3);

        if( (ptv.segment<3>(3) - G.block<3,1>(0, q_index)).norm() > 1e-5) {
          std::cout << __FILE__ << " " << __LINE__ << ": \n"
                    << ptv.transpose() << "\n"
                    << G.block<3,1>(0, q_index) << std::endl;
        }

      }
    } else {
      cerr << __FILE__ << " " << __LINE__ << ": "
           << "Custom joints not supported." << endl;
      abort();
      unsigned k = model.mJoints[j].custom_joint_index;
      G.block(0, q_index, 3, model.mCustomJoints[k]->mDoFCount) =
          ((point_trans
            * model.X_base[j].inverse()).toMatrix()
           * model.mCustomJoints[k]->S).block(
            3,0,3,model.mCustomJoints[k]->mDoFCount);
    }

    j = model.lambda[j];
  }
}

RBDL_DLLAPI void CalcPointJacobian6D (
    Model &model,
    ADModel &ad_model,
    VectorNd const &q,
    MatrixNd const &q_dirs,
    unsigned body_id,
    Vector3d const &point_position,
    MatrixNd &G,
    vector<MatrixNd> &ad_G,
    bool update_kinematics) {
  unsigned ndirs = q_dirs.cols();
  assert(ndirs == ad_G.size());

  if (update_kinematics) {
    UpdateKinematicsCustom (model, ad_model, &q, &q_dirs, 0, 0, 0, 0);
  }

  MatrixNd r_dirs(3, ndirs);
  Vector3d r = CalcBodyToBaseCoordinates(model, ad_model, q, q_dirs, body_id,
                                         point_position, r_dirs, false);

  vector<SpatialTransform> point_trans_dirs(ndirs);
  for (unsigned idir = 0; idir < ndirs; idir++) {
    point_trans_dirs[idir] = SpatialTransform(Matrix3d::Zero(), r_dirs.col(idir));
  }
  SpatialTransform point_trans = SpatialTransform (Matrix3d::Identity(), r);

  for (unsigned idir = 0; idir < ndirs; idir++) {
    assert(ad_G[idir].rows() == 6);
    assert(ad_G[idir].cols() == model.qdot_size);
  }
  assert (G.rows() == 6 && G.cols() == model.qdot_size );


  unsigned int reference_body_id = body_id;

  if (model.IsFixedBodyId(body_id)) {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
  }

  unsigned int j = reference_body_id;
  while (j != 0) {
    unsigned q_index = model.mJoints[j].q_index;

    if(model.mJoints[j].mJointType != JointTypeCustom){
      if (model.mJoints[j].mDoFCount == 1) {

        vector<SpatialTransform> Xbji_dirs(ndirs);
        SpatialTransform Xbji;
        inverse(ndirs,
                model.X_base[j], ad_model.X_base[j],
                Xbji, Xbji_dirs);

        vector<SpatialVector> XbjiSj_dirs(ndirs);
        SpatialVector XbjiSj;
        applySTSV(ndirs,
                  Xbji, Xbji_dirs,
                  model.S[j],
                  XbjiSj, XbjiSj_dirs);

        SpatialVector Gqi;
        vector<SpatialVector> Gqi_dirs(ndirs);
        applySTSV(ndirs,
                  point_trans, point_trans_dirs,
                  XbjiSj, XbjiSj_dirs,
                  Gqi, Gqi_dirs);

        for (unsigned idir = 0; idir < ndirs; idir++) {
          ad_G[idir].col(q_index) = Gqi_dirs[idir];
        }
        G.col(q_index) = Gqi;
      } else if (model.mJoints[j].mDoFCount == 3) {
        cerr << __FILE__ << " " << __LINE__
             << ": Multi-DoF joint not supported!" << endl;
        abort();
      }
    } else {
      cerr << __FILE__ << " " << __LINE__
           << ": Custom joint not supported!" << endl;
      abort();
    }
    j = model.lambda[j];
  }
}

// -----------------------------------------------------------------------------
} // namespace AD
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
