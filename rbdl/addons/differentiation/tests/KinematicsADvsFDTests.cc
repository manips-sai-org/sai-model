#include <UnitTest++.h>

#include <iostream>

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/rbdl_mathutils.h"

#include "KinematicsED.h"
#include "KinematicsAD.h"
#include "KinematicsFD.h"
#include "KinematicsFDC.h"

#include "Fixtures.h"
#include "Human36Fixture.h"

#include "ModelCheckADvsFD.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-8;

// -----------------------------------------------------------------------------

template<typename T>
void UpdateKinematicsCustomTemplate(
    T & obj,
    unsigned trial_count) {
  Model   ad_model = obj.model;
  Model   fd_model = obj.model;

  ADModel ad_d_model(ad_model);
  ADModel fd_d_model(fd_model);

  int const nq = obj.model.dof_count;
  unsigned ndirs = 3 * nq;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    VectorNd q    = VectorNd::Random(nq);
    VectorNd qd   = VectorNd::Random(nq);
    VectorNd qdd  = VectorNd::Random(nq);

    MatrixNd q_dirs = MatrixNd::Random(nq, ndirs);
    MatrixNd qd_dirs = MatrixNd::Random(nq, ndirs);
    MatrixNd qdd_dirs = MatrixNd::Random(nq, ndirs);

    RigidBodyDynamics::AD::UpdateKinematicsCustom(ad_model, ad_d_model, &q, &q_dirs, &qd, &qd_dirs, &qdd, &qdd_dirs);

    RigidBodyDynamics::FD::UpdateKinematicsCustom(fd_model, &fd_d_model, q, q_dirs, qd, qd_dirs, qdd, qdd_dirs);

    checkModelsADvsFD(
          ndirs,
          ad_model, ad_d_model,
          fd_model, fd_d_model);
  }
}

TEST_FIXTURE ( CartPendulum, CartPendulumUpdateKinematicsCustom) {
  UpdateKinematicsCustomTemplate(*this, 10);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXUpdateKinematicsCustom) {
  UpdateKinematicsCustomTemplate(*this, 10);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZUpdateKinematicsCustom) {
  UpdateKinematicsCustomTemplate(*this, 10);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpUpdateKinematicsCustom) {
  UpdateKinematicsCustomTemplate(*this, 10);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpUpdateKinematicsCustom) {
  UpdateKinematicsCustomTemplate(*this, 10);
}

TEST_FIXTURE ( Human36, Human36UpdateKinematicsCustom) {
  UpdateKinematicsCustomTemplate(*this, 10);
}

// -----------------------------------------------------------------------------

template<typename T>
void UpdateKinematicsCustomEDvsADTemplate(
    T & obj,
    unsigned const & trial_count
) {
  Model   ad_model = obj.model;
  Model   ed_model = obj.model;

  ADModel ad_d_model(ad_model);
  EDModel ed_d_model(ed_model);

  int const nq = obj.model.dof_count;
  unsigned ndirs = 3 * nq;
  ndirs = 1;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    VectorNd q    = VectorNd::Random(nq);
    VectorNd qd   = VectorNd::Random(nq);
    VectorNd qdd  = VectorNd::Random(nq);

    MatrixNd q_dirs = MatrixNd::Random(nq, ndirs);
    MatrixNd qd_dirs = MatrixNd::Random(nq, ndirs);
    MatrixNd qdd_dirs = MatrixNd::Random(nq, ndirs);

    RigidBodyDynamics::AD::UpdateKinematicsCustom(ad_model, ad_d_model, &q, &q_dirs, &qd, &qd_dirs, &qdd, &qdd_dirs);

    RigidBodyDynamics::ED::UpdateKinematicsCustom(ed_model, ed_d_model, &q, &q_dirs, &qd, &qd_dirs, &qdd, &qdd_dirs);

    checkModelsADvsED (
          ndirs,
          ad_model, ad_d_model,
          ed_model, ed_d_model);
  }
}

TEST_FIXTURE ( CartPendulum, CartPendulumUpdateKinematicsCustomEDvsAD ) {
  UpdateKinematicsCustomEDvsADTemplate(*this, 1);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXUpdateKinematicsCustomEDvsAD) {
  UpdateKinematicsCustomEDvsADTemplate(*this, 10);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZUpdateKinematicsCustomEDvsAD) {
  UpdateKinematicsCustomEDvsADTemplate(*this, 10);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpUpdateKinematicsCustomEDvsAD) {
  UpdateKinematicsCustomEDvsADTemplate(*this, 10);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpUpdateKinematicsCustomEDvsAD) {
  UpdateKinematicsCustomEDvsADTemplate(*this, 10);
}

TEST_FIXTURE ( Human36, Human36UpdateKinematicsCustomEDvsAD) {
  UpdateKinematicsCustomEDvsADTemplate(*this, 10);
}

// -----------------------------------------------------------------------------

template<typename T>
void UpdateKinematicsEDvsADTemplate(
    T & obj,
    unsigned const & trial_count
) {
  Model   ad_model = obj.model;
  Model   ed_model = obj.model;

  ADModel ad_d_model(ad_model);
  EDModel ed_d_model(ed_model);

  int const nq = obj.model.dof_count;
  unsigned ndirs = 3 * nq;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    VectorNd q    = VectorNd::Random(nq);
    VectorNd qd   = VectorNd::Random(nq);
    VectorNd qdd  = VectorNd::Random(nq);

    MatrixNd q_dirs = MatrixNd::Random(nq, ndirs);
    MatrixNd qd_dirs = MatrixNd::Random(nq, ndirs);
    MatrixNd qdd_dirs = MatrixNd::Random(nq, ndirs);

    RigidBodyDynamics::AD::UpdateKinematics(ad_model, ad_d_model, q, q_dirs, qd, qd_dirs, qdd, qdd_dirs);

    RigidBodyDynamics::ED::UpdateKinematics(ed_model, ed_d_model, q, q_dirs, qd, qd_dirs, qdd, qdd_dirs);

    checkModelsADvsED (
          ndirs,
          ad_model, ad_d_model,
          ed_model, ed_d_model);
  }
}

TEST_FIXTURE ( CartPendulum, CartPendulumUpdateKinematicsEDvsAD ) {
  UpdateKinematicsEDvsADTemplate(*this, 1);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXUpdateKinematicsEDvsAD) {
  UpdateKinematicsEDvsADTemplate(*this, 10);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZUpdateKinematicsEDvsAD) {
  UpdateKinematicsEDvsADTemplate(*this, 10);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpUpdateKinematicsEDvsAD) {
  UpdateKinematicsEDvsADTemplate(*this, 10);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpUpdateKinematicsEDvsAD) {
  UpdateKinematicsEDvsADTemplate(*this, 10);
}

TEST_FIXTURE ( Human36, Human36UpdateKinematicsEDvsAD) {
  UpdateKinematicsEDvsADTemplate(*this, 10);
}


// -----------------------------------------------------------------------------

template <typename T>
void CalcPointVelocityTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec) {
  Model   model = obj.model;
  ADModel ad_model(model);

  // set up input quantities
  int const nq = model.dof_count;
  unsigned const ndirs = 2 * nq;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    Vector3d point_position = Vector3d::Random();
    VectorNd q    = VectorNd::Random(nq);
    MatrixNd q_dirs = MatrixNd::Random(nq, ndirs);
    VectorNd qdot = VectorNd::Random(nq);
    MatrixNd qdot_dirs = MatrixNd::Random(nq, ndirs);

    // set up no output quantities
    Vector3d G    = Vector3d::Zero ();
    Vector3d ad_G = Vector3d::Zero ();
    Vector3d fd_G = Vector3d::Zero ();

    // set up derivative output quantities
    MatrixNd ad_G_dirs (3, ndirs);
    MatrixNd fd_G_dirs (3, ndirs);

    for (unsigned i = 1; i < model.mBodies.size(); i++) {
      unsigned int body_id = model.mBodyNameMap[model.GetBodyName(i)];
      if (body_id == 0) {
        continue;
      }
      // call nominal version
      G = CalcPointVelocity (
            model, q, qdot, body_id, point_position, true);

      // call FD version
      fd_G = RigidBodyDynamics::FD::CalcPointVelocity (
        model, q, q_dirs, qdot, qdot_dirs,
        body_id, point_position, fd_G_dirs
      );

      // call AD version
      ad_G = RigidBodyDynamics::AD::CalcPointVelocity (
        model, ad_model, q, q_dirs, qdot,
        qdot_dirs, body_id, point_position,
        ad_G_dirs, true
      );

      // nominal check
      CHECK_ARRAY_CLOSE(ad_G.data(), G.data(), G.size(), array_close_prec);
      CHECK_ARRAY_CLOSE(fd_G.data(), G.data(), G.size(), array_close_prec);
      CHECK_ARRAY_CLOSE(ad_G.data(), fd_G.data(), G.size(), array_close_prec);
      // derivative check
      CHECK_ARRAY_CLOSE(ad_G_dirs.data(), fd_G_dirs.data(),
                        3 * ndirs, array_close_prec);
    }
  }
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcPointVelocity) {
  CalcPointVelocityTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcPointVelocity) {
  CalcPointVelocityTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcPointVelocity) {
  CalcPointVelocityTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcPointVelocity) {
  CalcPointVelocityTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcPointVelocity) {
  CalcPointVelocityTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Human36, Human36CalcPointVelocity) {
  CalcPointVelocityTemplate(*this, 10, 1e-5);
}


// -----------------------------------------------------------------------------

template <typename T>
void CalcPointVelocity6DTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec) {
  Model   model = obj.model;
  Model   ad_model = obj.model;
  Model   fd_model = obj.model;
  ADModel ad_d_model(ad_model);
  ADModel fd_d_model(fd_model);

  // set up input quantities
  int const nq = ad_model.dof_count;
  unsigned const ndirs = 2 * nq;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    Vector3d pt_pos = Vector3d::Random();
    VectorNd q    = VectorNd::Random(nq);
    MatrixNd q_dirs = MatrixNd::Random(nq, ndirs);
    VectorNd qdot = VectorNd::Random(nq);
    MatrixNd qdot_dirs = MatrixNd::Random(nq, ndirs);

    // set up no output quantities
    SpatialVector pv6d;
    SpatialVector ad_pv6d;
    SpatialVector fd_pv6d;

    // set up derivative output quantities
    vector<SpatialVector> ad_pv6d_dirs (ndirs);
    vector<SpatialVector> fd_pv6d_dirs (ndirs);

    for (unsigned i = 1; i < model.mBodies.size(); i++) {
      unsigned int body_id = model.mBodyNameMap[ad_model.GetBodyName(i)];
      if (body_id == 0) {
        continue;
      }

      // call nominal version
      pv6d = CalcPointVelocity6D(model, q, qdot, body_id, pt_pos);

      // call FD version
      fd_pv6d = RigidBodyDynamics::FD::CalcPointVelocity6D (
            fd_model, &fd_d_model, q, q_dirs, qdot, qdot_dirs,
            body_id, pt_pos, fd_pv6d_dirs);

      // call AD version
      ad_pv6d = RigidBodyDynamics::AD::CalcPointVelocity6D (
            ad_model, ad_d_model, q, q_dirs, qdot, qdot_dirs,
            body_id, pt_pos, ad_pv6d_dirs);

      checkModelsADvsFD(ndirs,
                        ad_model, ad_d_model,
                        fd_model, fd_d_model);

      // nominal check
      CHECK_ARRAY_CLOSE(ad_pv6d.data(), pv6d.data(), 6,
                        array_close_prec);
      CHECK_ARRAY_CLOSE(fd_pv6d.data(), pv6d.data(), 6,
                        array_close_prec);
      // derivative check
      for (unsigned idir = 0; idir < ndirs; idir++) {
        CHECK_ARRAY_CLOSE(ad_pv6d_dirs[idir].data(), fd_pv6d_dirs[idir].data(),
                          6, array_close_prec);
      }
    }
  }
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcPointVelocity6D) {
  CalcPointVelocity6DTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcPointVelocity6D) {
  CalcPointVelocity6DTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcPointVelocity6D) {
  CalcPointVelocity6DTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcPointVelocity6D) {
  CalcPointVelocity6DTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcPointVelocity6D) {
  CalcPointVelocity6DTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Human36, Human36CalcPointVelocity6D) {
  CalcPointVelocity6DTemplate(*this, 10, 1e-5);
}


// -----------------------------------------------------------------------------

template <typename T>
void CalcBaseToBodyCoordinatesTemplate(
    T & obj,
    const unsigned trial_count,
    const double array_close_prec
) {
  Model   model = obj.model;
  Model   ad_model = obj.model;
  Model   fd_model = obj.model;
  ADModel ad_d_model(ad_model);
  ADModel fd_d_model(fd_model);

  // set up input quantities
  int const nq = model.dof_count;
  unsigned const ndirs = nq;
  bool update_kinematics = true;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    Vector3d point_position = Vector3d::Random();
    MatrixNd point_position_dirs = MatrixNd::Random(3, ndirs);

    VectorNd q = VectorNd::Random(nq);
    MatrixNd q_dirs = MatrixNd::Identity(nq, nq);

    // set up no output quantities
    Vector3d G    = Vector3d::Zero ();
    Vector3d ad_G = Vector3d::Zero ();
    Vector3d fd_G = Vector3d::Zero ();

    // set up derivative output quantities
    MatrixNd ad_G_dirs (3, ndirs);
    MatrixNd fd_G_dirs (3, ndirs);

    for (unsigned i = 1; i < model.mBodies.size(); i++) {
      unsigned int body_id = model.mBodyNameMap[model.GetBodyName(i)];
      if (body_id == 0) {
        continue;
      }

      G = CalcBaseToBodyCoordinates (
            model, q, body_id, point_position, update_kinematics);

      fd_G = RigidBodyDynamics::FD::CalcBaseToBodyCoordinates (
            fd_model, &fd_d_model, q, q_dirs, body_id, point_position,
            point_position_dirs, fd_G_dirs);

      ad_G = RigidBodyDynamics::AD::CalcBaseToBodyCoordinates (
            ad_model, ad_d_model, q, q_dirs, body_id, point_position,
            point_position_dirs, ad_G_dirs, update_kinematics);

      checkModelsADvsFD(ndirs,
                        ad_model, ad_d_model,
                        fd_model, fd_d_model);

      CHECK_ARRAY_CLOSE(ad_G.data(), G.data(), G.size(), array_close_prec);
      CHECK_ARRAY_CLOSE(fd_G.data(), G.data(), G.size(), array_close_prec);
      CHECK_ARRAY_CLOSE(ad_G_dirs.data(), fd_G_dirs.data(),
                        3 * ndirs, array_close_prec);
    }
  }
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcBaseToBodyCoordinates) {
  CalcBaseToBodyCoordinatesTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcBaseToBodyCoordinates) {
  CalcBaseToBodyCoordinatesTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcBaseToBodyCoordinates) {
  CalcBaseToBodyCoordinatesTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcBaseToBodyCoordinates) {
  CalcBaseToBodyCoordinatesTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcBaseToBodyCoordinates) {
  CalcBaseToBodyCoordinatesTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Human36, Human36CalcBaseToBodyCoordinates) {
  CalcBaseToBodyCoordinatesTemplate(*this, 10, 1e-5);
}

// -----------------------------------------------------------------------------

template <typename T>
void CalcBodyToBaseCoordinatesTemplate(
    T & obj,
    const unsigned trial_count,
    const double array_close_prec
) {
  Model ad_model = obj.model;
  Model fd_model = obj.model;
  ADModel ad_d_model(ad_model);
  ADModel fd_d_model(fd_model);

  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    bool update_kinematics = true;
    Vector3d point_position = Vector3d::Random();

    VectorNd q = VectorNd::Random(nq);
    MatrixNd q_dirs = MatrixNd::Identity(nq, nq);

    // set up no output quantities
    Vector3d G = Vector3d::Zero ();
    Vector3d G_ad = Vector3d::Zero ();
    Vector3d G_fd = Vector3d::Zero ();

    // set up derivative output quantities
    MatrixNd derivative_ad (3, ndirs);
    MatrixNd derivative_fd (3, ndirs);

    for (unsigned i = 1; i < ad_model.mBodies.size(); i++) {
      unsigned int body_id = ad_model.mBodyNameMap[ad_model.GetBodyName(i)];
      if (body_id == 0) {
        continue;
      }
      // call nominal version
      G = CalcBodyToBaseCoordinates (
            ad_model, q, body_id, point_position, update_kinematics
            );

      // call FD version
      G_fd = RigidBodyDynamics::FD::CalcBodyToBaseCoordinates (
            fd_model, &fd_d_model,
            q, q_dirs,
            body_id,
            point_position,
            derivative_fd
            );

      // call AD version
      G_ad = RigidBodyDynamics::AD::CalcBodyToBaseCoordinates (
            ad_model, ad_d_model,
            q, q_dirs,
            body_id,
            point_position,
            derivative_ad,
            update_kinematics
            );

      checkModelsADvsFD(ndirs,
                        ad_model, ad_d_model,
                        fd_model, fd_d_model);

      CHECK_ARRAY_CLOSE(G_ad.data(), G.data(), G.size(), array_close_prec);
      CHECK_ARRAY_CLOSE(G_fd.data(), G.data(), G.size(), array_close_prec);
      CHECK_ARRAY_CLOSE(G_ad.data(), G_fd.data(), G.size(), array_close_prec);

      CHECK_ARRAY_CLOSE(derivative_ad.data(), derivative_fd.data(),
                        3 * ndirs, array_close_prec);
    }
  }
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcBodyToBaseCoordinates) {
  CalcBodyToBaseCoordinatesTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcBodyToBaseCoordinates) {
  CalcBodyToBaseCoordinatesTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcBodyToBaseCoordinates) {
  CalcBodyToBaseCoordinatesTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcBodyToBaseCoordinates) {
  CalcBodyToBaseCoordinatesTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcBodyToBaseCoordinates) {
  CalcBodyToBaseCoordinatesTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Human36, Human36CalcBodyToBaseCoordinates) {
  CalcBodyToBaseCoordinatesTemplate(*this, 10, 1e-6);
}


// -----------------------------------------------------------------------------

template <typename T>
void CalcBodyWorldOrientationTemplate(
    T & obj,
    const unsigned trial_count,
    const double array_close_prec
) {
  Model & model = obj.model;
  ADModel ad_model(model);
  int const nq = model.dof_count;

  VectorNd q = VectorNd::Zero(nq);
  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);

  int ndirs = nq;
  vector<Matrix3d> ad_E_dirs(ndirs);
  vector<Matrix3d> fd_E_dirs(ndirs);

  unsigned trial = 0;
  do {
    for (unsigned i = 1; i < model.mBodies.size(); i++) {
      unsigned id = model.mBodyNameMap[model.GetBodyName(i)];
      if (id == 0) {
        continue;
      }

      Matrix3d ad_E = RigidBodyDynamics::AD::CalcBodyWorldOrientation(
        model, ad_model, q, q_dirs, id, ad_E_dirs
      );

      Matrix3d fd_E = RigidBodyDynamics::FD::CalcBodyWorldOrientation(
        model, q, q_dirs, id, fd_E_dirs
      );

      CHECK_ARRAY_CLOSE(ad_E.data(), fd_E.data(), 9, array_close_prec);
      for (int idir = 0; idir < ndirs; idir++) {
        CHECK_ARRAY_CLOSE(fd_E_dirs[idir].data(), ad_E_dirs[idir].data(), 9, array_close_prec);
      }
    }
    q = VectorNd::Random(nq);
  } while(trial++ < trial_count);
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcBodyWorldOrientation) {
  CalcBodyWorldOrientationTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcBodyWorldOrientation) {
  CalcBodyWorldOrientationTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcBodyWorldOrientation) {
  CalcBodyWorldOrientationTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcBodyWorldOrientation) {
  CalcBodyWorldOrientationTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcBodyWorldOrientation) {
  CalcBodyWorldOrientationTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Human36, Human36CalcBodyWorldOrientation) {
  CalcBodyWorldOrientationTemplate(*this, 10, 1e-6);
}


// -----------------------------------------------------------------------------

template <typename T>
void CalcPointAccelerationTemplate(
    T & obj,
    const unsigned trial_count,
    const double array_close_prec
) {
  Model   & model    = obj.model;
  ADModel ad_model(model);
  int const nq       = model.dof_count;

  VectorNd q    = VectorNd::Zero(nq);
  VectorNd qd   = VectorNd::Zero(nq);
  VectorNd qdd  = VectorNd::Zero(nq);
  MatrixNd q_dirs   = MatrixNd::Zero(nq, 3 * nq);
  MatrixNd qd_dirs  = MatrixNd::Zero(nq, 3 * nq);
  MatrixNd qdd_dirs = MatrixNd::Zero(nq, 3 * nq);

  q_dirs.block(0, 0, nq, nq)        = MatrixNd::Identity(nq, nq);
  qd_dirs.block(0, nq, nq, nq)      = MatrixNd::Identity(nq, nq);
  qdd_dirs.block(0, 2 * nq, nq, nq) = MatrixNd::Identity(nq, nq);

  int ndirs = 3 * nq;
  MatrixNd ad_a_dirs(3, ndirs);
  MatrixNd fd_a_dirs(3, ndirs);

  Vector3d reference_point = Vector3dZero;
  unsigned trial = 0;
  do {
    for (unsigned i = 1; i < model.mBodies.size(); i++) {
      unsigned id = model.mBodyNameMap[model.GetBodyName(i)];
      if (id == 0) {
        continue;
      }

      Vector3d ad_a = RigidBodyDynamics::AD::CalcPointAcceleration (
            model, ad_model,
            q, q_dirs, qd, qd_dirs, qdd, qdd_dirs,
            id, reference_point,
            ad_a_dirs);

      Vector3d fd_a = RigidBodyDynamics::FD::CalcPointAcceleration (
            model,
            q, q_dirs, qd, qd_dirs, qdd, qdd_dirs,
            id, reference_point, fd_a_dirs);

      CHECK_ARRAY_CLOSE(ad_a.data(), fd_a.data(), 3, array_close_prec);
      CHECK_ARRAY_CLOSE(fd_a_dirs.data(), ad_a_dirs.data(),
                        3 * ndirs, array_close_prec);
    }
    q.setRandom();
    qd.setRandom();
    qdd.setRandom();
    q_dirs.setRandom();
    qd_dirs.setRandom();
    qdd_dirs.setRandom();
  } while(trial++ < trial_count);
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcPointAcceleration) {
  CalcPointAccelerationTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcPointAcceleration) {
  CalcPointAccelerationTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcPointAcceleration) {
  CalcPointAccelerationTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcPointAcceleration) {
  CalcPointAccelerationTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcPointAcceleration) {
  CalcPointAccelerationTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE ( Human36, Human36CalcPointAcceleration) {
  CalcPointAccelerationTemplate(*this, 10, 1e-6);
}


// -----------------------------------------------------------------------------

template <typename T>
void CalcPointJacobianTemplate(
    T & obj,
    const unsigned trial_count,
    const double array_close_prec
) {
  Model model = obj.model;
  Model ad_model = obj.model;
  Model fd_model = obj.model;
  ADModel ad_d_model(ad_model);
  ADModel fd_d_model(fd_model);

  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;
  bool update_kinematics = true;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    Vector3d point_position = Vector3d::Random();

    VectorNd q = VectorNd::Random(nq);
    MatrixNd q_dirs = MatrixNd::Random(nq, nq);

    for (unsigned i = 2; i < ad_model.mBodies.size(); i++) {
      unsigned int body_id = ad_model.mBodyNameMap[ad_model.GetBodyName(i)];
      if (body_id == 0) {
        continue;
      }

      // set up no output quantities
      MatrixNd G = MatrixNd::Zero (3, nq);
      MatrixNd ad_G = MatrixNd::Zero (3, nq);
      MatrixNd fd_G = MatrixNd::Zero (3, nq);

      // set up derivative output quantities
      vector<MatrixNd> ad_G_dirs (ndirs, ad_G);
      vector<MatrixNd> fd_G_dirs (ndirs, fd_G);

      CalcPointJacobian (
        model, q, body_id, point_position, G, update_kinematics
      );

      RigidBodyDynamics::AD::CalcPointJacobian (
        ad_model, ad_d_model,
        q, q_dirs,
        body_id,
        point_position,
        ad_G, ad_G_dirs,
        update_kinematics
      );

      RigidBodyDynamics::FDC::CalcPointJacobian (
        fd_model, &fd_d_model,
        q, q_dirs,
        body_id,
        point_position,
        fd_G, fd_G_dirs
      );

      checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);

      CHECK_ARRAY_CLOSE(ad_G.data(), G.data(), G.size(), array_close_prec);
      CHECK_ARRAY_CLOSE(fd_G.data(), G.data(), G.size(), array_close_prec);

      for (unsigned idir = 0; idir < ndirs; idir++) {
        MatrixNd error;
        double max;
        error = (ad_G_dirs[idir] - fd_G_dirs[idir]).cwiseAbs();
        max = error.maxCoeff();
        if (max > array_close_prec) {
          std::cout << "error = \n" << error << std::endl;
          std::cout << "max   = " << max << std::endl;
          std::cout << "ad_G_dirs[" << idir << "] = \n"<< ad_G_dirs[idir] << std::endl;
          std::cout << "fd_G_dirs[" << idir << "] = \n"<< fd_G_dirs[idir] << std::endl;
          std::cout << endl;
        }

        CHECK_ARRAY_CLOSE(
              ad_G_dirs[idir].data(),
              fd_G_dirs[idir].data(),
              ad_G_dirs[idir].size(),
              array_close_prec);
      }
    }
  }
}


TEST_FIXTURE ( CartPendulum, CartPendulumCalcPointJacobian) {
  CalcPointJacobianTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcPointJacobian) {
  CalcPointJacobianTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcPointJacobian) {
  CalcPointJacobianTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcPointJacobian) {
  CalcPointJacobianTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcPointJacobian) {
  CalcPointJacobianTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE(
  MultiPendulumWithBranches,
  MultiPendulumWithBranchesCalcPointJacobian_4_2
) {
  this->create_model(4, 2);
  CalcPointJacobianTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE(
  MultiPendulumWithBranches,
  MultiPendulumWithBranchesCalcPointJacobian_6_3
) {
  this->create_model(6, 3);
  CalcPointJacobianTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE(
  MultiPendulumWithBranches,
  MultiPendulumWithBranchesCalcPointJacobian_12_3
) {
  this->create_model(12, 3);
  CalcPointJacobianTemplate(*this, 10, 1e-10);
}


TEST_FIXTURE ( Human36, Human36CalcPointJacobian) {
  CalcPointJacobianTemplate(*this, 10, 1e-9);
}

// -----------------------------------------------------------------------------

template <typename T>
void CalcPointJacobianEDvsADTemplate(
    T & obj,
    const unsigned trial_count,
    const double array_close_prec
) {
  Model model = obj.model;
  Model ad_model = obj.model;
  Model fd_model = obj.model;
  ADModel ad_d_model(ad_model);
  EDModel fd_d_model(fd_model);

  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;
  bool update_kinematics = true;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    Vector3d point_position = Vector3d::Random();

    VectorNd q = VectorNd::Random(nq);
    MatrixNd q_dirs = MatrixNd::Random(nq, nq);

    for (unsigned i = 2; i < ad_model.mBodies.size(); i++) {
      unsigned int body_id = ad_model.mBodyNameMap[ad_model.GetBodyName(i)];
      if (body_id == 0) {
        continue;
      }

      // set up no output quantities
      MatrixNd G = MatrixNd::Zero (3, nq);
      MatrixNd ad_G = MatrixNd::Zero (3, nq);
      MatrixNd fd_G = MatrixNd::Zero (3, nq);

      // set up derivative output quantities
      vector<MatrixNd> ad_G_dirs (ndirs, ad_G);
      vector<MatrixNd> fd_G_dirs (ndirs, fd_G);

      CalcPointJacobian (
        model, q, body_id, point_position, G, update_kinematics
      );

      RigidBodyDynamics::AD::CalcPointJacobian (
        ad_model, ad_d_model,
        q, q_dirs,
        body_id,
        point_position,
        ad_G, ad_G_dirs,
        update_kinematics
      );

      RigidBodyDynamics::ED::CalcPointJacobian (
        fd_model, fd_d_model,
        q, q_dirs,
        body_id,
        point_position,
        fd_G, fd_G_dirs
      );

      checkModelsADvsED(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);

      CHECK_ARRAY_CLOSE(ad_G.data(), G.data(), G.size(), array_close_prec);
      CHECK_ARRAY_CLOSE(fd_G.data(), G.data(), G.size(), array_close_prec);

      for (unsigned idir = 0; idir < ndirs; idir++) {
        MatrixNd error;
        double max;
        error = (ad_G_dirs[idir] - fd_G_dirs[idir]).cwiseAbs();
        max = error.maxCoeff();
        if (max > array_close_prec) {
          std::cout << "error = \n" << error << std::endl;
          std::cout << "max   = " << max << std::endl;
          std::cout << "ad_G_dirs[" << idir << "] = \n"<< ad_G_dirs[idir] << std::endl;
          std::cout << "ed_G_dirs[" << idir << "] = \n"<< fd_G_dirs[idir] << std::endl;
          std::cout << endl;
        }

        CHECK_ARRAY_CLOSE(
          ad_G_dirs[idir].data(),
          fd_G_dirs[idir].data(),
          ad_G_dirs[idir].size(),
          array_close_prec
        );
      }
    }
  }
}

TEST_FIXTURE (CartPendulum, CartPendulumCalcPointJacobianEDvsAD)
{
  CalcPointJacobianEDvsADTemplate (*this, 1, 1e-10);
}

TEST_FIXTURE (Arm2DofX, Arm2DofXCalcPointJacobianEDvsAD)
{
  CalcPointJacobianEDvsADTemplate (*this, 10, 1e-10);
}


TEST_FIXTURE (Arm2DofZ, Arm2DofZCalcPointJacobianEDvsAD)
{
  CalcPointJacobianEDvsADTemplate (*this, 10, 1e-10);
}

TEST_FIXTURE (Arm3DofXZYp, Arm3DofXZYpCalcPointJacobianEDvsAD)
{
  CalcPointJacobianEDvsADTemplate (*this, 10, 1e-10);
}

TEST_FIXTURE (Arm3DofXZZp, Arm3DofXZZpCalcPointJacobianEDvsAD)
{
  CalcPointJacobianEDvsADTemplate (*this, 10, 1e-10);
}

TEST_FIXTURE (
  MultiPendulumWithBranches,
  MultiPendulumWithBranchesCalcPointJacobianEDvsAD_4_2
) {
  this->create_model(4, 2);
  CalcPointJacobianEDvsADTemplate (*this, 10, 1e-10);
}

TEST_FIXTURE (
  MultiPendulumWithBranches,
  MultiPendulumWithBranchesCalcPointJacobianEDvsAD_6_3
) {
  this->create_model(6, 3);
  CalcPointJacobianEDvsADTemplate (*this, 10, 1e-10);
}

TEST_FIXTURE (
  MultiPendulumWithBranches,
  MultiPendulumWithBranchesCalcPointJacobianEDvsAD_12_3
) {
  this->create_model(12, 3);
  CalcPointJacobianEDvsADTemplate (*this, 10, 1e-10);
}

TEST_FIXTURE (Human36, Human36CalcPointJacobianEDvsAD)
{
  CalcPointJacobianEDvsADTemplate (*this, 10, 1e-9);
}


// -----------------------------------------------------------------------------
template <typename T>
void CalcPointJacobian6DTemplate(
    T & obj,
    const unsigned trial_count,
    const double array_close_prec
) {
  Model model = obj.model;
  Model ad_model = obj.model;
  Model fd_model = obj.model;
  ADModel ad_d_model(ad_model);
  ADModel fd_d_model(fd_model);

  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;
  bool update_kinematics = true;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    Vector3d point_position = Vector3d::Random();

    VectorNd q = VectorNd::Random(nq);
    MatrixNd q_dirs = MatrixNd::Random(nq, nq);

    for (unsigned i = 2; i < ad_model.mBodies.size(); i++) {
      unsigned int body_id = ad_model.mBodyNameMap[ad_model.GetBodyName(i)];
      if (body_id == 0) {
        continue;
      }

      // set up no output quantities
      MatrixNd G = MatrixNd::Zero (6, nq);
      MatrixNd ad_G = MatrixNd::Zero (6, nq);
      MatrixNd fd_G = MatrixNd::Zero (6, nq);

      // set up derivative output quantities
      vector<MatrixNd> ad_G_dirs (ndirs, ad_G);
      vector<MatrixNd> fd_G_dirs (ndirs, fd_G);

      CalcPointJacobian6D (
            model, q, body_id, point_position, G, update_kinematics);

      RigidBodyDynamics::AD::CalcPointJacobian6D (
            ad_model, ad_d_model,
            q, q_dirs,
            body_id,
            point_position,
            ad_G, ad_G_dirs,
            update_kinematics);

      RigidBodyDynamics::FDC::CalcPointJacobian6D (
            fd_model, &fd_d_model,
            q, q_dirs,
            body_id,
            point_position,
            fd_G, fd_G_dirs);

      checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);

      CHECK_ARRAY_CLOSE(ad_G.data(), G.data(), G.rows() * G.cols(), array_close_prec);
      CHECK_ARRAY_CLOSE(fd_G.data(), G.data(), G.rows() * G.cols(), array_close_prec);

      for (unsigned idir = 0; idir < ndirs; idir++) {
        MatrixNd error;
        double max;
        error = (ad_G_dirs[idir] - fd_G_dirs[idir]).cwiseAbs();
        max = error.maxCoeff();
        if (max > array_close_prec) {
          std::cout << "error = \n" << error << std::endl;
          std::cout << "max   = " << max << std::endl;
          std::cout << "ad_G_dirs[" << idir << "] = \n"<< ad_G_dirs[idir] << std::endl;
          std::cout << "fd_G_dirs[" << idir << "] = \n"<< fd_G_dirs[idir] << std::endl;
          std::cout << endl;
        }

        CHECK_ARRAY_CLOSE(
              ad_G_dirs[idir].data(),
              fd_G_dirs[idir].data(),
              ad_G_dirs[idir].rows() * ad_G_dirs[idir].cols(),
              array_close_prec);
      }
    }
  }
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcPointJacobian6D) {
  CalcPointJacobian6DTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcPointJacobian6D) {
  CalcPointJacobian6DTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcPointJacobian6D) {
  CalcPointJacobian6DTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcPointJacobian6D) {
  CalcPointJacobian6DTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcPointJacobian6D) {
  CalcPointJacobian6DTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE(
  MultiPendulumWithBranches,
  MultiPendulumWithBranchesCalcPointJacobian6D_4_2
) {
  this->create_model(4, 2);
  CalcPointJacobian6DTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE(
  MultiPendulumWithBranches,
  MultiPendulumWithBranchesCalcPointJacobian6D_6_3
) {
  this->create_model(6, 3);
  CalcPointJacobian6DTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE(
  MultiPendulumWithBranches,
  MultiPendulumWithBranchesCalcPointJacobian6D_12_3
) {
  this->create_model(12, 3);
  CalcPointJacobian6DTemplate(*this, 10, 1e-10);
}

TEST_FIXTURE ( Human36, Human36CalcPointJacobian6D) {
  CalcPointJacobian6DTemplate(*this, 10, 1e-9);
}


// -----------------------------------------------------------------------------

