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

#include "JointED.h"
#include "ModelED.h"
#include "KinematicsED.h"
#include "SpatialAlgebraOperatorsED.h"

#include "rbdl_mathutilsAD.h"
#include "rbdl_mathutilsFD.h"

using std::cerr;
using std::endl;
using std::vector;

using namespace RigidBodyDynamics::Math;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace ED {
// -----------------------------------------------------------------------------

/*
RBDL_DLLAPI Vector3d CalcBodyToBaseCoordinates (
    Model &model,
    EDModel &ad_model,
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
    ED::UpdateKinematicsCustom (
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
    EDModel & ad_model,
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
*/

RBDL_DLLAPI void UpdateKinematicsCustom (
    Model &model,
    EDModel & ed_model,
    const VectorNd * q,
    const MatrixNd * q_dirs,
    const VectorNd * qdot,
    const MatrixNd * qdot_dirs,
    const VectorNd * qddot,
    const MatrixNd * qddot_dirs
) {
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
    assert(q_dirs && q_dirs->cols() == ndirs);
  }

  ed_model.resize_directions(ndirs);

  if (q) {
    for (i = 1; i < model.mBodies.size(); i++) {
      unsigned int q_index = model.mJoints[i].q_index;
      unsigned int lambda = model.lambda[i];

      VectorNd QDot_zero (VectorNd::Zero (model.q_size));

      jcalc (model, i, (*q), QDot_zero);

      // nominal evaluation
      // model.X_lambda[i] = model.X_J[i] * model.X_T[i];  // removed
      // derivative evaluation
      for (unsigned int idir = 0; idir < ndirs; ++idir)
      {
        SpatialVector v = model.S[i]*q_dirs->row(q_index)[idir];
        ed_model.X_lambda[i][idir] = SpatialTransform(
          -VectorCrossMatrix(v.head(3)) * model.X_lambda[i].E,
          v.tail(3).transpose() * model.X_lambda[i].E
        );
      }

      if (lambda != 0) {
        // derivative evaluation
        for (unsigned int idir = 0; idir < ndirs; ++idir)
        {
          ed_model.X_base[i][idir] = SpatialTransform (
                // E * XT.E,
                ed_model.X_lambda[i][idir].E * model.X_base[lambda].E
                + model.X_lambda[i].E * ed_model.X_base[lambda][idir].E,
                // XT.r + XT.E.transpose() * r
                ed_model.X_base[lambda][idir].r
                  + ed_model.X_base[lambda][idir].E.transpose() * model.X_lambda[i].r
                  + model.X_base[lambda].E.transpose() * ed_model.X_lambda[i][idir].r
          );
        }
        // nominal evaluation
        model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
      } else {
        // derivative evaluation
        for (unsigned int idir = 0; idir < ndirs; ++idir)
        {
          ed_model.X_base[i][idir] = ed_model.X_lambda[i][idir];
        }
        // nominal evaluation
        model.X_base[i] = model.X_lambda[i];
      }
    }
  }

  if (qdot) {
    for (i = 1; i < model.mBodies.size(); i++) {
      unsigned int q_index = model.mJoints[i].q_index;
      unsigned int lambda = model.lambda[i];

      jcalc (model, i, *q, *qdot);

      if (lambda != 0) {
        // nominal evaluation
        model.v[i] = model.X_lambda[i].apply(model.v[lambda]);

        // derivative evaluation
        ed_model.v[i].leftCols(ndirs)
            = crossm(model.v[i])*model.S[i]*q_dirs->row(q_index)
            + model.X_lambda[i].toMatrix()*ed_model.v[lambda].leftCols(ndirs)
            + model.S[i]*qdot_dirs->row(q_index);
        // nominal evaluation continued
        model.v[i] += model.v_J[i];

        // nominal evaluation
        model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);
        // derivative evaluation
        ed_model.c[i].leftCols(ndirs) =
            crossm(model.v[i])*model.S[i]*qdot_dirs->row(q_index)
            - crossm(model.v_J[i]) * (ed_model.v[i].leftCols(ndirs) );

      } else {
        // derivative evaluation
        ed_model.v[i].leftCols(ndirs) = model.S[i]*qdot_dirs->row(q_index);
        // nominal evaluation
        model.v[i] = model.v_J[i];

        // nominal evaluation
        model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);
        // derivative evaluation
        ed_model.c[i].leftCols(ndirs) =
            crossm(model.v[i])*model.S[i]*qdot_dirs->row(q_index)
            - crossm(model.v_J[i]) * (ed_model.v[i].leftCols(ndirs) );
      }
    }
  }

  if (qddot) {
    for (i = 1; i < model.mBodies.size(); i++) {
      unsigned int q_index = model.mJoints[i].q_index;
      unsigned int lambda = model.lambda[i];

      if (lambda != 0) {
        SpatialVector a = model.X_lambda[i].apply(model.a[lambda]);

        // derivative evaluation
        ed_model.a[i].leftCols(ndirs)
        = crossm(a)*model.S[i]*(*q_dirs).row(model.mJoints[i].q_index) +
          model.X_lambda[i].toMatrix()*ed_model.a[lambda].leftCols(ndirs)
          + ed_model.c[i].leftCols(ndirs);

        // nominal evaluation
        model.a[i] = a + model.c[i];
      } else {
        // derivative evaluation
        ed_model.a[i].leftCols(ndirs) = ed_model.c[i].leftCols(ndirs);
        // nominal evaluation
        model.a[i] = model.c[i];
      }

      if( model.mJoints[i].mJointType != JointTypeCustom){
        if (model.mJoints[i].mDoFCount == 1) {
          // derivative evaluation
          ed_model.a[i].leftCols(ndirs) += model.S[i] * (*qddot_dirs).row(q_index);
          // nominal evaluation
          model.a[i] = model.a[i] + model.S[i] * (*qddot)[q_index];
        } else if (model.mJoints[i].mDoFCount == 3) {
          std::cerr << "NOT SUPPORTED" << std::endl;
          abort();
        }
      } else {
          std::cerr << "NOT SUPPORTED" << std::endl;
          abort();
      }
    }
  }
}

RBDL_DLLAPI void UpdateKinematics (
    Model &model,
    EDModel &ed_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs,
    const VectorNd &qddot,
    const MatrixNd &qddot_dirs
) {
  LOG << "-------- " << __func__ << " --------" << std::endl;

  unsigned int i = 0;
  const unsigned int ndirs = q_dirs.cols();
  ed_model.resize_directions(ndirs);

  // derivative evaluation
  for (unsigned int idir = 0; idir < ndirs; ++idir) {
    ed_model.a[0].setZero();
  }
  // nominal evaluation
  model.a[0].setZero();

  for (i = 1; i < model.mBodies.size(); i++) {
    unsigned int q_index = model.mJoints[i].q_index;
    unsigned int lambda = model.lambda[i];

    // nominal evaluation
    // NOTE a derivative is analytically computed later on
    jcalc (model, i, q, qdot);

    // nominal evaluation
    // model.X_lambda[i] = model.X_J[i] * model.X_T[i];  // removed
    // derivative evaluation
    for (unsigned int idir = 0; idir < ndirs; ++idir)
    {
      SpatialVector v = model.S[i]*q_dirs.row(q_index)[idir];
      ed_model.X_lambda[i][idir] = SpatialTransform(
        -VectorCrossMatrix(v.head(3)) * model.X_lambda[i].E,
        v.tail(3).transpose() * model.X_lambda[i].E
      );
    }

    if (lambda != 0) {
      // derivative evaluation
      for (unsigned int idir = 0; idir < ndirs; ++idir)
      {
        ed_model.X_base[i][idir] = SpatialTransform (
              // E * XT.E,
              ed_model.X_lambda[i][idir].E * model.X_base[lambda].E
              + model.X_lambda[i].E * ed_model.X_base[lambda][idir].E,
              // XT.r + XT.E.transpose() * r
              ed_model.X_base[lambda][idir].r
                + ed_model.X_base[lambda][idir].E.transpose() * model.X_lambda[i].r
                + model.X_base[lambda].E.transpose() * ed_model.X_lambda[i][idir].r
        );
      }
      // nominal evaluation
      model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];

      // nominal evaluation
      model.v[i] = model.X_lambda[i].apply(model.v[lambda]);
      // derivative evaluation
      ed_model.v[i].leftCols(ndirs)
          = crossm(model.v[i])*model.S[i]*q_dirs.row(q_index)
          + model.X_lambda[i].toMatrix()*ed_model.v[lambda].leftCols(ndirs)
          + model.S[i]*qdot_dirs.row(q_index);

      // nominal evaluation continued
      model.v[i] += model.v_J[i];

    } else {
      // derivative evaluation
      for (unsigned int idir = 0; idir < ndirs; ++idir)
      {
        ed_model.X_base[i][idir] = ed_model.X_lambda[i][idir];
      }
      // nominal evaluation
      model.X_base[i] = model.X_lambda[i];

      // derivative evaluation
      ed_model.v[i].leftCols(ndirs) = model.S[i]*qdot_dirs.row(q_index);
      // nominal evaluation
      model.v[i] = model.v_J[i];
    }

    // nominal evaluation
    model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);
    // derivative evaluation
    ed_model.c[i].leftCols(ndirs) =
        crossm(model.v[i])*model.S[i]*qdot_dirs.row(model.mJoints[i].q_index)
        - crossm(model.v_J[i]) * (ed_model.v[i].leftCols(ndirs) );

    // nominal evaluation
    // model.a[i] = model.X_lambda[i].apply(model.a[lambda]) + model.c[i];
    model.a[i] = model.X_lambda[i].apply(model.a[lambda]);
    // derivative evaluation
    ed_model.a[i] = crossm(model.a[i])*model.S[i]*q_dirs.row(q_index)
        + model.X_lambda[i].toMatrix()*ed_model.a[lambda].leftCols(ndirs)
        + ed_model.c[i].leftCols(ndirs);
    // nominal evaluation continued
    model.a[i] += model.c[i];

    if(model.mJoints[i].mJointType != JointTypeCustom){
      if (model.mJoints[i].mDoFCount == 1) {
        // derivative evaluation
        ed_model.a[i] += model.S[i] * qddot_dirs.row(q_index);
        // nominal evaluation
        model.a[i] = model.a[i] + model.S[i] * qddot[q_index];
      } else if (model.mJoints[i].mDoFCount == 3) {
        cerr << "Multi-dof not supported." << endl;
        abort();
      }
    } else {
      cerr << __FILE__ << " " << __LINE__
           << ": Custom joints not supported." << endl;
      abort();
    }
  }
}


RBDL_DLLAPI Matrix3d CalcBodyWorldOrientation (
    Model & model,
    EDModel & ed_model,
    VectorNd const & q,
    MatrixNd const & q_dirs,
    const unsigned int body_id,
    vector<Matrix3d> & ad_derivative,
    bool update_kinematics
) {
  const unsigned ndirs = q_dirs.cols();
  assert(ad_derivative.size() == ndirs);

  // update the Kinematics if necessary
  if (update_kinematics) {
    RigidBodyDynamics::ED::UpdateKinematicsCustom (model, ed_model, &q, &q_dirs, 0, 0, 0, 0);
  }

  if (body_id >= model.fixed_body_discriminator) {
    std::cerr << "Fixed bodies not yet supported!" << std::endl;
    abort();
  }

  for (unsigned int idir = 0; idir < ndirs; idir++) {
    ad_derivative[idir] = ed_model.X_base[body_id][idir].E;
  }
  // nominal evaluation
  return model.X_base[body_id].E;
}


/*
RBDL_DLLAPI Vector3d CalcPointVelocity (
    Model & model,
    EDModel & ad_model,
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
    EDModel &ad_model,
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
*/

RBDL_DLLAPI Vector3d CalcPointAcceleration (
    Model &model,
    EDModel &ed_model,
    const Math::VectorNd &q,
    const Math::MatrixNd &q_dirs,
    const Math::VectorNd &qdot,
    const Math::MatrixNd &qdot_dirs,
    const Math::VectorNd &qddot,
    const Math::MatrixNd &qddot_dirs,
    unsigned int body_id,
    const Math::Vector3d &point_position,
    Math::MatrixNd &ed_derivative,
    bool update_kinematics
) {
  LOG << "-------- " << __func__ << " --------" << std::endl;

  int const ndirs = q_dirs.cols();
//  assert(ndirs == ed_derivative.cols());
  assert(3     == ed_derivative.rows());

  ed_model.resize_directions(ndirs);

  // Reset the velocity of the root body
  model.v[0].setZero();
  model.a[0].setZero();

  if (update_kinematics) {
    // derivative evaluation
    RigidBodyDynamics::ED::UpdateKinematics (
      model, ed_model, q, q_dirs, qdot, qdot_dirs, qddot, qddot_dirs
    );
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

  Matrix3d E = RigidBodyDynamics::ED::CalcBodyWorldOrientation(
    model, ed_model, q, q_dirs, reference_body_id, ed_model.ad_E, false
  );

  // derivative evaluation
  // TODO directly build up SpatialTransform
  for (int idir = 0; idir < ndirs; idir++) {
    Matrix3d ad_ETrx = ed_model.ad_E[idir].transpose() * Matrix3d(
          0., -reference_point[2],  reference_point[1],
        reference_point[2],                  0., -reference_point[0],
        -reference_point[1],  reference_point[0],                  0.);
    // + E.transpose() * Matrix3dZero;
    ed_model.ad_p_X_i[idir].block<3,3>(0, 0) = ed_model.ad_E[idir].transpose();
    ed_model.ad_p_X_i[idir].block<3,3>(3, 0) = -ad_ETrx;
    ed_model.ad_p_X_i[idir].block<3,3>(3, 3) = ed_model.ad_E[idir].transpose();
  }
  // nominal evaluation
  SpatialTransform p_X_i (E.transpose(), reference_point);

  // derivative evaluation
  for (int idir = 0; idir < ndirs; idir++) {
    ed_model.ad_p_v_i.col(idir) = ed_model.ad_p_X_i[idir] * model.v[reference_body_id]
        + p_X_i.apply(ed_model.v[reference_body_id].col(idir));
  }
  // nominal evaluation
  SpatialVector p_v_i = p_X_i.apply(model.v[reference_body_id]);

  // derivative evaluation
  ed_model.ad_a_dash.leftCols(ndirs) =
      RigidBodyDynamics::AD::cross3d(p_v_i.segment<3>(0))*ed_model.ad_p_v_i.bottomRows(3)
      - RigidBodyDynamics::AD::cross3d(p_v_i.segment<3>(3))*ed_model.ad_p_v_i.topRows(3);
  // nominal evaluation
  Vector3d a_dash = p_v_i.segment<3>(0).cross(p_v_i.segment<3>(3));

  // derivative evaluation
  for (int idir = 0; idir < ndirs; idir++) {
    ed_derivative.col(idir)
        = (ed_model.ad_p_X_i[idir] * model.a[reference_body_id] + p_X_i.apply(ed_model.a[reference_body_id].col(idir))).bottomRows(3);
  }
  ed_derivative.leftCols(ndirs) += ed_model.ad_a_dash;
  // nominal evaluation
  return p_X_i.apply(model.a[reference_body_id]).segment<3>(3) + a_dash;
}

RBDL_DLLAPI
void CalcPointJacobian (
    Model &model,
    EDModel &ed_model,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    unsigned body_id,
    const Vector3d &point_position,
    MatrixNd &G,
    vector<MatrixNd> &G_dirs,
    bool update_kinematics
) {
  LOG << "-------- " << __func__ << " --------" << std::endl;

  const unsigned int ndirs = q_dirs.cols();
  ed_model.resize_directions(ndirs);

  // update the Kinematics if necessary
  if (update_kinematics)
  {
    UpdateKinematicsCustom (
      model, ed_model,
      &q, &q_dirs,
      nullptr, nullptr,
      nullptr, nullptr
    );
  }

  // derivative evaluation
  for (unsigned idir = 0; idir < ndirs; idir++) {
    // NOTE point_body_coordinates is constant, no derivative needed!
    ed_model.X_0[idir] = SpatialTransform (
        Matrix3d::Zero(),
        // code from CalcBodyToBaseCoordinates
        ed_model.X_base[body_id][idir].r
        + ed_model.X_base[body_id][idir].E.transpose() * point_position
    );
  }

  // nominal evaluation
  SpatialTransform point_trans = SpatialTransform (
      Matrix3d::Identity(),
      CalcBodyToBaseCoordinates (
        model,
        q,
        body_id,
        point_position,
        false
      )
    );

  // derivative evaluation
  for (unsigned int idirs = 0; idirs < ndirs; idirs++) {
    assert (G_dirs[idirs].rows() == 3);
    assert (G_dirs[idirs].cols() == model.qdot_size );
  }
  // nominal evaluation
  assert (G.rows() == 3 && G.cols() == model.qdot_size );

  unsigned int reference_body_id = body_id;

  if (model.IsFixedBodyId(body_id))
  {
    unsigned int fbody_id = body_id - model.fixed_body_discriminator;
    reference_body_id = model.mFixedBodies[fbody_id].mMovableParent;
  }

  unsigned int j = reference_body_id;

  if (body_id >= model.fixed_body_discriminator) {
    cerr << __FILE__ << " " << __LINE__
         << ": Fixed bodies not yet supported!" << endl;
    abort();
  }


  // e[j] is set to 1 if joint j contributes to the Jacobian that we are
  // computing. For all other joints the column will be zero.
  while (j != 0)
  {
    unsigned int q_index = model.mJoints[j].q_index;

    if(model.mJoints[j].mJointType != JointTypeCustom)
    {
      if (model.mJoints[j].mDoFCount == 1)
      {
        // nominal evaluation
        SpatialVector v = model.X_base[j].inverse().apply(model.S[j]);
        SpatialVector w = point_trans.apply(v);
        SpatialVector dv, dw;
        // derivative evaluation
        for (unsigned idir = 0; idir < ndirs; idir++)
        {
          // derivative evaluation
          dv.head(3) = ed_model.X_base[j][idir].E.transpose() * model.S[j].head(3);
          dv.tail(3) = ed_model.X_base[j][idir].E.transpose() * (
                model.S[j].tail<3>() - model.S[j].head<3>().cross(model.X_base[j].E * model.X_base[j].r)
            ) + model.X_base[j].E.transpose() * (
              ed_model.X_base[j][idir].E * model.X_base[j].r
              + model.X_base[j].E * ed_model.X_base[j][idir].r
            ).cross(model.S[j].head<3>());

          // derivative evaluation
          dw.head(3) = ed_model.X_0[idir].E * v.head(3) + point_trans.E * dv.head(3);
          dw.tail(3) = point_trans.E * (
              dv.tail<3>()
              + dv.head<3>().cross(point_trans.r)
              + v.head<3>().cross(ed_model.X_0[idir].r)
            ) - ed_model.X_0[idir].E * point_trans.r.cross(v.head<3>());

          // derivative evaluation
          G_dirs[idir].block<3, 1>(0, q_index) = dw.tail(3);
        }
        // nominal evaluation
        G.block(0,q_index, 3, 1) = w.tail<3>();
      }
      else if (model.mJoints[j].mDoFCount == 3)
      {
        cerr << __FILE__ << " " << __LINE__ << ": "
             << "Multi-DoF joints not supported." << endl;
        abort();
      }
    }
    else
    {
      cerr << __FILE__ << " " << __LINE__ << ": "
           << "Custom joints not supported." << endl;
      abort();
    }

    j = model.lambda[j];
  }
};

/*
RBDL_DLLAPI void CalcPointJacobian6D (
    Model &model,
    EDModel &ed_model,
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
*/

// -----------------------------------------------------------------------------
} // namespace ED
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
