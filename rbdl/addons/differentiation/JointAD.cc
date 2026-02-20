/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <limits>
#include <assert.h>

#include "rbdl/Logging.h"

#include "rbdl/Model.h"
#include "rbdl/Joint.h"

#include "rbdl_mathutilsAD.h"

#include "JointAD.h"
#include "ModelAD.h"

#include "SpatialAlgebraOperatorsAD.h"


using std::cerr;
using std::endl;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

using namespace RigidBodyDynamics::Math;

RBDL_DLLAPI void jcalc (
    Model & model,
    ADModel & ad_model,
    unsigned int joint_id,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs
) {
  unsigned const ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  ad_model.resize_directions(ndirs);

  // exception if we calculate it for the root body
  assert (joint_id > 0);

  if (model.mJoints[joint_id].mJointType == JointTypeRevoluteX) {
    // derivative code
    for (unsigned idir = 0; idir < ndirs; ++idir) {
      /// TODO: Most entries of q_dirs are ZERO !!!
      /// => shorten that!

      AD::Xrotx (
            q[model.mJoints[joint_id].q_index],
            q_dirs(model.mJoints[joint_id].q_index, idir),
            ad_model.X_J[joint_id][idir]);

      //            ad_model.X_J[joint_id][idir] = Math::AD::Xrotx(
      //								q[model.mJoints[joint_id].q_index],
      //								q_dirs(model.mJoints[joint_id].q_index, idir));
      ad_model.v_J[joint_id][idir][0] = qdot_dirs(model.mJoints[joint_id].q_index, idir); // v_J = S*qdot
      ad_model.c_J[joint_id][idir].setZero(); // c_J = Sdot*qdot
    }

    // nominal code
    model.X_J[joint_id] = Math::Xrotx (q[model.mJoints[joint_id].q_index]);
    model.v_J[joint_id][0] = qdot[model.mJoints[joint_id].q_index];
  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteY) {
    // derivative code
    for (unsigned idir = 0; idir < ndirs; ++idir) {
      /// TODO: Most entries of q_dirs are ZERO !!!
      /// => shorten that!
      AD::Xroty(
            q[model.mJoints[joint_id].q_index],
            q_dirs(model.mJoints[joint_id].q_index, idir),
            ad_model.X_J[joint_id][idir]);
      //            ad_model.X_J[joint_id][idir] = Math::AD::Xroty(
      //								q[model.mJoints[joint_id].q_index],
      //								q_dirs(model.mJoints[joint_id].q_index, idir));
      ad_model.v_J[joint_id][idir][1] = qdot_dirs(model.mJoints[joint_id].q_index, idir); // v_J = S*qdot
      ad_model.c_J[joint_id][idir].setZero(); // c_J = Sdot*qdot
    }
    // nominal code
    model.X_J[joint_id] = Math::Xroty (q[model.mJoints[joint_id].q_index]);
    model.v_J[joint_id][1] = qdot[model.mJoints[joint_id].q_index];
  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteZ) {
    // derivative code
    for (unsigned int idir = 0; idir < ndirs; ++idir) {
      /// TODO: Most entries of q_dirs are ZERO !!!
      /// => shorten that!
      AD::Xrotz (
            q[model.mJoints[joint_id].q_index],
          q_dirs(model.mJoints[joint_id].q_index, idir),
          ad_model.X_J[joint_id][idir]);
      //            ad_model.X_J[joint_id][idir] = Math::AD::Xrotz(
      //								q[model.mJoints[joint_id].q_index],
      //								q_dirs(model.mJoints[joint_id].q_index, idir));
      ad_model.v_J[joint_id][idir][2] = qdot_dirs(model.mJoints[joint_id].q_index, idir); // v_J = S*qdot
      ad_model.c_J[joint_id][idir].setZero(); // c_J = Sdot*qdot
    }
    // nominal code
    model.X_J[joint_id] = Math::Xrotz (q[model.mJoints[joint_id].q_index]);
    model.v_J[joint_id][2] = qdot[model.mJoints[joint_id].q_index];
  } else if (model.S[joint_id] == SpatialVector (0., 0., 0., 1., 0., 0.)) {
    // derivative code
    for (unsigned int idir = 0; idir < ndirs; ++idir) {
      ad_model.X_J[joint_id][idir] = AD::Xtrans(
            Vector3d (q(model.mJoints[joint_id].q_index), 0.0, 0.0),
            Vector3d (q_dirs(model.mJoints[joint_id].q_index, idir), 0.0, 0.0)
            );
      ad_model.v_J[joint_id][idir][3] = qdot_dirs(model.mJoints[joint_id].q_index, idir); // v_J = S*qdot
      ad_model.c_J[joint_id][idir] = SpatialVector::Zero(); // v_J = Sdot*qdot
    }
    // nominal code
    model.X_J[joint_id] = Math::Xtrans (Vector3d (q[model.mJoints[joint_id].q_index], 0., 0.));
    model.v_J[joint_id][3] = qdot[model.mJoints[joint_id].q_index];
  } else if (model.mJoints[joint_id].mDoFCount == 1) {
    // S = [0,. 0., 0., a, b, c] a,b,c in {0,1}, a + b + c = 1
    int vidx = -1;
    for (int i = 0; i < 6; i++) {
      if (model.S[joint_id][i] == 1.) {
        vidx = i;
      }
    }
    for (unsigned int idir = 0; idir < ndirs; ++idir) {
      // NOTE: X_T is a constant model dependent transformation
      ad_model.X_J[joint_id][idir] = jcalc_XJ (model, ad_model, joint_id, idir, q, q_dirs);
      ad_model.v_J[joint_id][idir][vidx] = qdot_dirs(model.mJoints[joint_id].q_index, idir); // v_J = S*qdot
      ad_model.c_J[joint_id][idir].setZero(); // v_J = Sdot*qdot
    }
    // nominal code
    model.X_J[joint_id] = jcalc_XJ (model, joint_id, q);
    model.v_J[joint_id] = model.S[joint_id] * qdot[model.mJoints[joint_id].q_index];
  } else if (model.mJoints[joint_id].mJointType == JointTypeSpherical) {
    std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerZYX) {
    std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerXYZ) {
    std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerYXZ) {
    std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeTranslationXYZ) {
    std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
    abort();
  } else {
    std::cerr << "Error: invalid joint type " << model.mJoints[joint_id].mJointType << " at id " << joint_id << std::endl;
    abort();
  }

  mulSTST(ndirs,
          model.X_J[joint_id], ad_model.X_J[joint_id],
          model.X_T[joint_id],
          model.X_lambda[joint_id], ad_model.X_lambda[joint_id]);
}


RBDL_DLLAPI Math::SpatialTransform jcalc_XJ (
    Model & model,
    ADModel & ad_model,
    unsigned int joint_id,
    unsigned int idir,
    const Math::VectorNd & q,
    const Math::MatrixNd & q_dirs) {
  // exception if we calculate it for the root body
  assert (joint_id > 0);

  if (model.mJoints[joint_id].mDoFCount == 1) {
    if (model.mJoints[joint_id].mJointType == JointTypeRevolute) {
      std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
      abort();
      // TODO derive Xrot in ad_Xrot
      // return Xrot (q[model.mJoints[joint_id].q_index], Vector3d (
      // 	model.mJoints[joint_id].mJointAxes[0][0],
      // 	model.mJoints[joint_id].mJointAxes[0][1],
      // 	model.mJoints[joint_id].mJointAxes[0][2]
      // 	));
    } else if (model.mJoints[joint_id].mJointType == JointTypePrismatic) {
      return AD::Xtrans(
            Vector3d (
              model.mJoints[joint_id].mJointAxes[0][3] * q(model.mJoints[joint_id].q_index),
          model.mJoints[joint_id].mJointAxes[0][4] * q(model.mJoints[joint_id].q_index),
          model.mJoints[joint_id].mJointAxes[0][5] * q(model.mJoints[joint_id].q_index)),
          Vector3d (
          model.mJoints[joint_id].mJointAxes[0][3] * q_dirs(model.mJoints[joint_id].q_index, idir),
          model.mJoints[joint_id].mJointAxes[0][4] * q_dirs(model.mJoints[joint_id].q_index, idir),
          model.mJoints[joint_id].mJointAxes[0][5] * q_dirs(model.mJoints[joint_id].q_index, idir)));
    }
  }
  std::cerr << "Error: invalid joint type!" << std::endl;
  abort();
  return SpatialTransform();
}

RBDL_DLLAPI void jcalc_X_lambda_S (
    Model &model,
    ADModel &ad_model,
    unsigned int joint_id,
    const VectorNd &q,
    const MatrixNd &q_dirs
    ) {
  // exception if we calculate it for the root body
  assert (joint_id > 0);

  size_t ndirs = q_dirs.cols();
  ad_model.resize_directions(ndirs);

  if (model.mJoints[joint_id].mJointType == JointTypeRevoluteX) {
    // derivative evaluation
//    for (unsigned int idir = 0; idir < ndirs; ++idir) {
//      // NOTE: X_T is a constant model dependent transformation
//      ad_model.X_lambda[joint_id][idir] = AD::Xrotx (
//            q[model.mJoints[joint_id].q_index],
//          q_dirs(model.mJoints[joint_id].q_index, idir))
//          * model.X_T[joint_id];
//    }
    SpatialTransform st_xrotx = Math::Xrotx (q[model.mJoints[joint_id].q_index]);
    std::vector<SpatialTransform> ad_st_xrotx(ndirs);
    for (unsigned idir = 0; idir < ndirs; ++idir) {
      ad_st_xrotx[idir] = AD::Xrotx (
          q[model.mJoints[joint_id].q_index],
          q_dirs(model.mJoints[joint_id].q_index, idir));
    }
    mulSTST(ndirs, st_xrotx, ad_st_xrotx,
        model.X_T[joint_id],
        model.X_lambda[joint_id],
        ad_model.X_lambda[joint_id]);
    // nominal evaluation
    model.X_lambda[joint_id] = st_xrotx * model.X_T[joint_id];
    model.S[joint_id] = model.mJoints[joint_id].mJointAxes[0];
  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteY) {
//    // derivative evaluation
//    for (unsigned int idir = 0; idir < ndirs; ++idir) {
//      // NOTE: X_T is a constant model dependent transformation
//      ad_model.X_lambda[joint_id][idir] = AD::Xroty (
//            q[model.mJoints[joint_id].q_index],
//          q_dirs(model.mJoints[joint_id].q_index, idir))
//          * model.X_T[joint_id];
//    }
    SpatialTransform st_xroty = Math::Xroty (q[model.mJoints[joint_id].q_index]);
    std::vector<SpatialTransform> ad_st_xroty(ndirs);
    for (unsigned idir = 0; idir < ndirs; ++idir) {
        ad_st_xroty[idir] = AD::Xroty (
            q[model.mJoints[joint_id].q_index],
            q_dirs(model.mJoints[joint_id].q_index, idir));
    }
    mulSTST(ndirs, st_xroty, ad_st_xroty,
          model.X_T[joint_id],
          model.X_lambda[joint_id],
          ad_model.X_lambda[joint_id]);
    // nominal evaluation
    model.X_lambda[joint_id] = st_xroty * model.X_T[joint_id];
    model.S[joint_id] = model.mJoints[joint_id].mJointAxes[0];
  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteZ) {
    // derivative evaluation
//    for (unsigned idir = 0; idir < ndirs; ++idir) {
//      // NOTE: X_T is a constant model dependent transformation
//      ad_model.X_lambda[joint_id][idir] = AD::Xrotz (
//            q[model.mJoints[joint_id].q_index],
//          q_dirs(model.mJoints[joint_id].q_index, idir))
//          * model.X_T[joint_id];
//    }
    SpatialTransform st_xrotz = Math::Xrotz (q[model.mJoints[joint_id].q_index]);
    std::vector<SpatialTransform> ad_st_xrotz(ndirs);
    for (unsigned idir = 0; idir < ndirs; ++idir) {
        ad_st_xrotz[idir] = AD::Xrotz (
                    q[model.mJoints[joint_id].q_index],
                  q_dirs(model.mJoints[joint_id].q_index, idir));
    }
    mulSTST(ndirs, st_xrotz, ad_st_xrotz,
            model.X_T[joint_id],
            model.X_lambda[joint_id],
            ad_model.X_lambda[joint_id]);
    // nominal evaluation
    model.X_lambda[joint_id] = st_xrotz * model.X_T[joint_id];
    model.S[joint_id] = model.mJoints[joint_id].mJointAxes[0];
  } else if (model.mJoints[joint_id].mDoFCount == 1) {
    // derivative evaluation
    for (unsigned int idir = 0; idir < ndirs; ++idir) {
      // NOTE: X_T is a constant model dependent transformation
      ad_model.X_lambda[joint_id][idir] = jcalc_XJ (model, ad_model,
                                                    joint_id, idir, q, q_dirs) * model.X_T[joint_id];
    }
    // nominal evaluation
    model.X_lambda[joint_id] = jcalc_XJ (model, joint_id, q) * model.X_T[joint_id];
    model.S[joint_id] = model.mJoints[joint_id].mJointAxes[0];
  } else if (model.mJoints[joint_id].mJointType == JointTypeSpherical) {
    std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerZYX) {
    std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerXYZ) {
    std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerYXZ ) {
    std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeTranslationXYZ ) {
    std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
    abort();
  } else {
    std::cerr << "Error: invalid joint type!" << std::endl;
    abort();
  }
}

// -----------------------------------------------------------------------------
} /* AD */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------
