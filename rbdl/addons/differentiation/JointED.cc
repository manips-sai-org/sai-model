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

#include "JointED.h"
#include "ModelED.h"
#include "ModelAD.h"
#include "SpatialAlgebraOperatorsED.h"

#include <cmath>

#ifdef __APPLE__
inline void sincosp(RigidBodyDynamics::Math::Scalar x, RigidBodyDynamics::Math::Scalar* sinp, RigidBodyDynamics::Math::Scalar* cosp) {
    return __sincos(x, sinp, cosp);
}
#else

inline void sincosp(RigidBodyDynamics::Math::Scalar x, RigidBodyDynamics::Math::Scalar* sinp, RigidBodyDynamics::Math::Scalar* cosp) {
    *sinp = sin(x);
    *cosp = cos(x);
}
#endif

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace ED {
// -----------------------------------------------------------------------------

using std::cerr;
using std::endl;
using namespace RigidBodyDynamics::Math;

RBDL_DLLAPI void jcalc (
    RigidBodyDynamics::Model &model,
    RigidBodyDynamics::EDModel &ed_model,
    unsigned int joint_id,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs
) {
  // exception if we calculate it for the root body
  assert (joint_id > 0);

  // receive number of directions
  const unsigned ndirs = q_dirs.cols();
  if(ndirs != qdot_dirs.cols()){abort();};
  // if(ndirs != qddot_dirs.cols()){abort()};

  // resize container if necessary
  ed_model.resize_directions(ndirs);

if (model.mJoints[joint_id].mJointType == JointTypeRevoluteX) {
    Scalar s, c;
    sincosp (q[model.mJoints[joint_id].q_index], &s, &c);

    model.X_lambda[joint_id].E = Matrix3d (
        model.X_T[joint_id].E(0, 0),
        model.X_T[joint_id].E(0, 1),
        model.X_T[joint_id].E(0, 2),

         c * model.X_T[joint_id].E(1, 0) + s * model.X_T[joint_id].E(2, 0),
         c * model.X_T[joint_id].E(1, 1) + s * model.X_T[joint_id].E(2, 1),
         c * model.X_T[joint_id].E(1, 2) + s * model.X_T[joint_id].E(2, 2),

        -s * model.X_T[joint_id].E(1, 0) + c * model.X_T[joint_id].E(2, 0),
        -s * model.X_T[joint_id].E(1, 1) + c * model.X_T[joint_id].E(2, 1),
        -s * model.X_T[joint_id].E(1, 2) + c * model.X_T[joint_id].E(2, 2));

    model.X_lambda[joint_id].r = model.X_T[joint_id].r;

    model.v_J[joint_id][0] = qdot[model.mJoints[joint_id].q_index];
  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteY) {
    Scalar s, c;
    sincosp (q[model.mJoints[joint_id].q_index], &s, &c);

    model.X_lambda[joint_id].E = Matrix3d (
        c * model.X_T[joint_id].E(0, 0) + -s * model.X_T[joint_id].E(2, 0),
        c * model.X_T[joint_id].E(0, 1) + -s * model.X_T[joint_id].E(2, 1),
        c * model.X_T[joint_id].E(0, 2) + -s * model.X_T[joint_id].E(2, 2),

        model.X_T[joint_id].E(1, 0),
        model.X_T[joint_id].E(1, 1),
        model.X_T[joint_id].E(1, 2),

        s * model.X_T[joint_id].E(0, 0) + c * model.X_T[joint_id].E(2, 0),
        s * model.X_T[joint_id].E(0, 1) + c * model.X_T[joint_id].E(2, 1),
        s * model.X_T[joint_id].E(0, 2) + c * model.X_T[joint_id].E(2, 2));

    model.X_lambda[joint_id].r = model.X_T[joint_id].r;

    model.v_J[joint_id][1] = qdot[model.mJoints[joint_id].q_index];
  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteZ) {
    Scalar s, c;
    sincosp (q[model.mJoints[joint_id].q_index], &s, &c);

    model.X_lambda[joint_id].E = Matrix3d (
         c * model.X_T[joint_id].E(0, 0) + s * model.X_T[joint_id].E(1, 0),
         c * model.X_T[joint_id].E(0, 1) + s * model.X_T[joint_id].E(1, 1),
         c * model.X_T[joint_id].E(0, 2) + s * model.X_T[joint_id].E(1, 2),

        -s * model.X_T[joint_id].E(0, 0) + c * model.X_T[joint_id].E(1, 0),
        -s * model.X_T[joint_id].E(0, 1) + c * model.X_T[joint_id].E(1, 1),
        -s * model.X_T[joint_id].E(0, 2) + c * model.X_T[joint_id].E(1, 2),

        model.X_T[joint_id].E(2, 0),
        model.X_T[joint_id].E(2, 1),
        model.X_T[joint_id].E(2, 2));

    model.X_lambda[joint_id].r = model.X_T[joint_id].r;

    model.v_J[joint_id][2] = qdot[model.mJoints[joint_id].q_index];
  } else if (model.mJoints[joint_id].mDoFCount == 1 &&
      model.mJoints[joint_id].mJointType != JointTypeCustom) {
      SpatialTransform X_J = jcalc_XJ (model, joint_id, q);
          model.v_J[joint_id] =
            model.S[joint_id] * qdot[model.mJoints[joint_id].q_index];
          model.X_lambda[joint_id] = X_J * model.X_T[joint_id];
  } else if (model.mJoints[joint_id].mJointType == JointTypeSpherical) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerZYX) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerXYZ) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerYXZ) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if(model.mJoints[joint_id].mJointType == JointTypeTranslationXYZ){
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeCustom) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Custom joint not supported." << endl;
    abort();
  } else {
    std::cerr << "Error: invalid joint type " << model.mJoints[joint_id].mJointType << " at id " << joint_id << std::endl;
    abort();
  }

  // model.X_lambda[joint_id] = model.X_J[joint_id] * model.X_T[joint_id];
  return;
}

// RBDL_DLLAPI void jcalc (
//     RigidBodyDynamics::Model &model,
//     RigidBodyDynamics::EDModel &ed_model,
//     unsigned int joint_id,
//     const VectorNd &q,
//     const MatrixNd &q_dirs,
//     const VectorNd &qdot,
//     const MatrixNd &qdot_dirs
// ) {
//   // exception if we calculate it for the root body
//   assert (joint_id > 0);

//   // receive number of directions
//   const unsigned ndirs = q_dirs.cols();
//   if(ndirs != qdot_dirs.cols()){abort();};
//   // if(ndirs != qddot_dirs.cols()){abort()};

//   // resize container if necessary
//   ed_model.resize_directions(ndirs);

//   if (model.mJoints[joint_id].mJointType == JointTypeRevoluteX) {
//     model.X_J[joint_id] = Xrotx (q[model.mJoints[joint_id].q_index]);
//     model.v_J[joint_id][0] = qdot[model.mJoints[joint_id].q_index];
//   } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteY) {
//     model.X_J[joint_id] = Xroty (q[model.mJoints[joint_id].q_index]);
//     model.v_J[joint_id][1] = qdot[model.mJoints[joint_id].q_index];
//   } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteZ) {
//     model.X_J[joint_id] = Xrotz (q[model.mJoints[joint_id].q_index]);
//     model.v_J[joint_id][2] = qdot[model.mJoints[joint_id].q_index];
//   } else if (model.mJoints[joint_id].mDoFCount == 1 &&
//       model.mJoints[joint_id].mJointType != JointTypeCustom) {
//     model.X_J[joint_id] = jcalc_XJ (model, joint_id, q);
//     model.v_J[joint_id]
//       = model.S[joint_id] * qdot[model.mJoints[joint_id].q_index];
//   } else if (model.mJoints[joint_id].mJointType == JointTypeSpherical) {
//     cerr << __FILE__ << " " << __LINE__ << ":"
//          << " Multi-DoF joint not supported." << endl;
//     abort();
//   } else if (model.mJoints[joint_id].mJointType == JointTypeEulerZYX) {
//     cerr << __FILE__ << " " << __LINE__ << ":"
//          << " Multi-DoF joint not supported." << endl;
//     abort();
//   } else if (model.mJoints[joint_id].mJointType == JointTypeEulerXYZ) {
//     cerr << __FILE__ << " " << __LINE__ << ":"
//          << " Multi-DoF joint not supported." << endl;
//     abort();
//   } else if (model.mJoints[joint_id].mJointType == JointTypeEulerYXZ) {
//     cerr << __FILE__ << " " << __LINE__ << ":"
//          << " Multi-DoF joint not supported." << endl;
//     abort();
//   } else if(model.mJoints[joint_id].mJointType == JointTypeTranslationXYZ){
//     cerr << __FILE__ << " " << __LINE__ << ":"
//          << " Multi-DoF joint not supported." << endl;
//     abort();
//   } else if (model.mJoints[joint_id].mJointType == JointTypeCustom) {
//     cerr << __FILE__ << " " << __LINE__ << ":"
//          << " Custom joint not supported." << endl;
//     abort();
//   } else {
//     std::cerr << "Error: invalid joint type " << model.mJoints[joint_id].mJointType << " at id " << joint_id << std::endl;
//     abort();
//   }

//   model.X_lambda[joint_id] = model.X_J[joint_id] * model.X_T[joint_id];
//   return;
// }

/*
RBDL_DLLAPI void jcalc (
    RigidBodyDynamics::Model &model,
    RigidBodyDynamics::EDModel &ed_model,
    unsigned int joint_id,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs
) {

  // exception if we calculate it for the root body
  assert (joint_id > 0);

  // receive number of directions
  const unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  // assert(ndirs == qddot_dirs.cols());

  // resize container if necessary
  ed_model.resize_directions(ndirs);

  if (model.mJoints[joint_id].mJointType == JointTypeRevoluteX) {
    // derivative evaluation
    Math::ED::Xrotx(
      model.X_J[joint_id], ed_model.X_J[joint_id],
      q[model.mJoints[joint_id].q_index],
      q_dirs.row(model.mJoints[joint_id].q_index),
      ndirs
    );
    // nominal evaluation
    // NOTE nominal evaluation happens in derivative code to save cos/sin operations
    // model.X_J[joint_id] = Xrotx (q[model.mJoints[joint_id].q_index]);

    // derivative evaluation
    ed_model.v_J[joint_id].leftCols(ndirs).row(0) = qdot_dirs.row(model.mJoints[joint_id].q_index);
    // nominal evaluation
    model.v_J[joint_id][0] = qdot[model.mJoints[joint_id].q_index];

  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteY) {

    // derivative evaluation
    Math::ED::Xroty(
      model.X_J[joint_id], ed_model.X_J[joint_id],
      q[model.mJoints[joint_id].q_index],
      q_dirs.row(model.mJoints[joint_id].q_index),
      ndirs
    );
    // nominal evaluation
    // NOTE nominal evaluation happens in derivative code to save cos/sin operations
    // model.X_J[joint_id] = Xroty (q[model.mJoints[joint_id].q_index]);

    // derivative evaluation
    ed_model.v_J[joint_id].leftCols(ndirs).row(1) = qdot_dirs.row(model.mJoints[joint_id].q_index);
    // nominal evaluation
    model.v_J[joint_id][1] = qdot[model.mJoints[joint_id].q_index];

  } else if (model.mJoints[joint_id].mJointType == JointTypeRevoluteZ) {

    // derivative evaluation
    // NOTE we can save common evaluations of cos and sin
    Math::ED::Xrotz(
      model.X_J[joint_id], ed_model.X_J[joint_id],
      q[model.mJoints[joint_id].q_index],
      q_dirs.row(model.mJoints[joint_id].q_index),
      ndirs
    );
    // // nominal evaluation
    // NOTE nominal evaluation happens in derivative code to save cos/sin operations
    // model.X_J[joint_id] = Xrotz (q[model.mJoints[joint_id].q_index]);

    // derivative evaluation
    ed_model.v_J[joint_id].leftCols(ndirs).row(2) = qdot_dirs.row(model.mJoints[joint_id].q_index);
    // nominal evaluation
    model.v_J[joint_id][2] = qdot[model.mJoints[joint_id].q_index];

  } else if (model.mJoints[joint_id].mDoFCount == 1 &&
      model.mJoints[joint_id].mJointType != JointTypeCustom) {
    // S = [0,. 0., 0., a, b, c] a,b,c in {0,1}, a + b + c = 1
    // derivative evaluation
    RigidBodyDynamics::ED::jcalc_XJ (
      model, ed_model,
      joint_id, ndirs,
      q, q_dirs,
      model.X_J[joint_id],
      ed_model.X_J[joint_id]
    );
    // nominal evaluation
    // model.X_J[joint_id] = jcalc_XJ (model, joint_id, q);

    // derivative evaluation
    int vidx = -1;
    for (int i = 0; i < 6; i++) {
      if (model.S[joint_id][i] == 1.) {
        vidx = i;
      }
    }
    ed_model.v_J[joint_id].leftCols(ndirs).row(vidx) = qdot_dirs.row(model.mJoints[joint_id].q_index);
    // nominal evaluation
    model.v_J[joint_id][vidx] = qdot[model.mJoints[joint_id].q_index];

  } else if (model.mJoints[joint_id].mJointType == JointTypeSpherical) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerZYX) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerXYZ) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeEulerYXZ) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if(model.mJoints[joint_id].mJointType == JointTypeTranslationXYZ){
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Multi-DoF joint not supported." << endl;
    abort();
  } else if (model.mJoints[joint_id].mJointType == JointTypeCustom) {
    cerr << __FILE__ << " " << __LINE__ << ":"
         << " Custom joint not supported." << endl;
    abort();
  } else {
    std::cerr << "Error: invalid joint type " << model.mJoints[joint_id].mJointType << " at id " << joint_id << std::endl;
    abort();
  }

  // derivative evaluation
  for (unsigned idir = 0; idir < ndirs; idir++) {
    ed_model.X_lambda[joint_id][idir] = ed_model.X_J[joint_id][idir] * model.X_T[joint_id];
  }
  // nominal evaluation
  model.X_lambda[joint_id] = model.X_J[joint_id] * model.X_T[joint_id];

  return;
}
*/

// RBDL_DLLAPI
// void jcalc_XJ (
//   RigidBodyDynamics::Model &model,
//   RigidBodyDynamics::EDModel &ed_model,
//   const unsigned int &joint_id,
//   const unsigned int &ndir,
//   const Math::VectorNd &q,
//   const Math::MatrixNd &q_dirs,
//   SpatialTransform &X,
//   std::vector<SpatialTransform> &X_dir
// ) {
//   assert (joint_id > 0);

//   if (model.mJoints[joint_id].mDoFCount == 1
//       && model.mJoints[joint_id].mJointType != JointTypeCustom) {
//     if (model.mJoints[joint_id].mJointType == JointTypeRevolute) {
//       return Xrot (q[model.mJoints[joint_id].q_index], Vector3d (
//             model.mJoints[joint_id].mJointAxes[0][0],
//             model.mJoints[joint_id].mJointAxes[0][1],
//             model.mJoints[joint_id].mJointAxes[0][2]
//             ));
//     } else if (model.mJoints[joint_id].mJointType == JointTypePrismatic) {
//       return Xtrans ( Vector3d (
//             model.mJoints[joint_id].mJointAxes[0][3]
//             * q[model.mJoints[joint_id].q_index],
//             model.mJoints[joint_id].mJointAxes[0][4]
//             * q[model.mJoints[joint_id].q_index],
//             model.mJoints[joint_id].mJointAxes[0][5]
//             * q[model.mJoints[joint_id].q_index]
//             )
//           );
//     }
//   }
//   std::cerr << "Error: invalid joint type!" << std::endl;
//   abort();
//   return SpatialTransform();
// }

/*
RBDL_DLLAPI
void jcalc_XJ (
  RigidBodyDynamics::Model &model,
  RigidBodyDynamics::EDModel &ed_model,
  const unsigned int &joint_id,
  const unsigned int &ndir,
  const Math::VectorNd &q,
  const Math::MatrixNd &q_dirs,
  SpatialTransform &X,
  std::vector<SpatialTransform> &X_dir
) {
  // exception if we calculate it for the root body
  assert (joint_id > 0);

  if (model.mJoints[joint_id].mDoFCount == 1
      && model.mJoints[joint_id].mJointType != JointTypeCustom) {
    if (model.mJoints[joint_id].mJointType == JointTypeRevolute) {
      std::cerr << "In: " << __func__ << "this joint type is not supported." << std::endl;
      abort();
      return;
    } else if (model.mJoints[joint_id].mJointType == JointTypePrismatic) {
      // derivative evaluation
      for (unsigned int idir = 0; idir < ndir; ++idir) {
        X_dir[idir].E.setZero ();
        X_dir[idir].r = Vector3d (
            model.mJoints[joint_id].mJointAxes[0][3]
            * q_dirs(model.mJoints[joint_id].q_index, idir),
            model.mJoints[joint_id].mJointAxes[0][4]
            * q_dirs(model.mJoints[joint_id].q_index, idir),
            model.mJoints[joint_id].mJointAxes[0][5]
            * q_dirs(model.mJoints[joint_id].q_index, idir)
        );
      }
      // nominal evaluation
      X.E = Matrix3d::Identity(3, 3);
      X.r = Vector3d (
        model.mJoints[joint_id].mJointAxes[0][3]
        * q[model.mJoints[joint_id].q_index],
        model.mJoints[joint_id].mJointAxes[0][4]
        * q[model.mJoints[joint_id].q_index],
        model.mJoints[joint_id].mJointAxes[0][5]
        * q[model.mJoints[joint_id].q_index]
      );
    }
  }
}
*/

// -----------------------------------------------------------------------------
} /* ED */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------
