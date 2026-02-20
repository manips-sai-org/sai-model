#include <UnitTest++.h>

#include <iostream>

#include "rbdl/Model.h"
#include "rbdl/Dynamics.h"
#include "rbdl/rbdl_mathutils.h"

#include "Fixtures.h"
#include "Human36Fixture.h"
#include "ModelAD.h"
#include "ModelED.h"
#include "DynamicsAD.h"
#include "DynamicsED.h"
#include "DynamicsFD.h"
#include "DynamicsFDC.h"
#include "rbdl_mathutilsAD.h"

#include "ModelCheckADvsFD.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-8;

// -----------------------------------------------------------------------------


template <typename T>
void ForwardDynamicsADTestTemplate(
    T & obj,
    unsigned int trial_count,
    double array_close_prec
) {
  Model ad_model = obj.model;
  Model fd_model = obj.model;
  Model fdc_model = obj.model;

  ADModel ad_d_model(ad_model);
  ADModel fd_d_model(fd_model);
  ADModel fdc_d_model(fdc_model);

  srand(666);

  for(unsigned int trial = 0; trial < trial_count; trial++) {
    VectorNd q = VectorNd::Random(obj.model.q_size);
    VectorNd qdot = VectorNd::Random(obj.model.q_size);
    VectorNd tau = VectorNd::Random(obj.model.q_size);

    unsigned int ndirs = 3 * obj.model.q_size;
    MatrixNd x = MatrixNd::Identity(ndirs, ndirs);
    MatrixNd q_dirs = x.block(0, 0, obj.model.q_size, ndirs);
    MatrixNd qdot_dirs = x.block(obj.model.q_size, 0, obj.model.q_size, ndirs);
    MatrixNd tau_dirs = x.block(2 * obj.model.q_size, 0, obj.model.q_size, ndirs);

    vector<SpatialVector> f_ext (
          obj.model.mBodies.size(),
          SpatialVector::Zero()
          );

    vector<vector<SpatialVector>> f_ext_dirs (
          obj.model.mBodies.size(),
          vector<SpatialVector>(
            ndirs,
            SpatialVector::Zero()));

    if (trial != 0) {
      for (unsigned i = 0; i < f_ext.size(); i++) {
        f_ext[i].setRandom();
      }
      for (unsigned i = 0; i < f_ext.size(); i++) {
        for (unsigned idir = 0; idir < f_ext_dirs[i].size(); idir++) {
          f_ext_dirs[i][idir].setRandom();
        }
      }
    }

    VectorNd ad_qddot (VectorNd::Zero(obj.model.q_size));
    VectorNd fd_qddot (VectorNd::Zero(obj.model.q_size));
    VectorNd fdc_qddot (VectorNd::Zero(obj.model.q_size));
    MatrixNd ad_dqddot  = MatrixNd::Zero(obj.model.qdot_size, ndirs);
    MatrixNd fd_dqddot  = MatrixNd::Zero(obj.model.qdot_size, ndirs);
    MatrixNd fdc_dqddot = MatrixNd::Zero(obj.model.qdot_size, ndirs);

    AD::ForwardDynamics(
          ad_model,
          ad_d_model,
          q, q_dirs,
          qdot, qdot_dirs,
          tau, tau_dirs,
          ad_qddot, ad_dqddot,
          &f_ext, &f_ext_dirs);

    FD::ForwardDynamics(
          fd_model,
          &fd_d_model,
          q, q_dirs,
          qdot, qdot_dirs,
          tau, tau_dirs,
          fd_qddot, fd_dqddot,
          &f_ext, &f_ext_dirs);

    FDC::ForwardDynamics(
          fdc_model,
          &fdc_d_model,
          q, q_dirs,
          qdot, qdot_dirs,
          tau, tau_dirs,
          fdc_qddot, fdc_dqddot,
          &f_ext, &f_ext_dirs
    );

    checkModelsADvsFD(
          ndirs,
          ad_model, ad_d_model,
          fdc_model, fdc_d_model);

    CHECK_ARRAY_CLOSE(
      ad_qddot.data(),
      fdc_qddot.data(),
      obj.model.q_size,
      array_close_prec
    );

    CHECK_ARRAY_CLOSE(
      fdc_dqddot.data(),
      ad_dqddot.data(),
      fd_dqddot.cols() * fd_dqddot.rows(),
      array_close_prec
    );

    checkModelsADvsFD(
          ndirs,
          ad_model, ad_d_model,
          fd_model, fd_d_model);

    CHECK_ARRAY_CLOSE(
      ad_qddot.data(),
      fd_qddot.data(),
      obj.model.q_size,
      array_close_prec
    );

    CHECK_ARRAY_CLOSE(
      fd_dqddot.data(),
      ad_dqddot.data(),
      fd_dqddot.cols() * fd_dqddot.rows(),
      array_close_prec
    );
  }
}

TEST_FIXTURE(CartPendulum, CartPendulumForwardDynamicsADTest){
  ForwardDynamicsADTestTemplate(*this, 1, 1e-5);
}

TEST_FIXTURE(Arm2DofX, Arm2DofXForwardDynamicsADTest){
  ForwardDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm2DofZ, Arm2DofZForwardDynamicsADTest){
  ForwardDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm3DofXZYp, Arm3DofXZYpForwardDynamicsADTest){
  ForwardDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm3DofXZZp, Arm3DofXZZpForwardDynamicsADTest){
  ForwardDynamicsADTestTemplate(*this, 10, 1e-5);
}

//TEST_FIXTURE( Human36, Human36ForwardDynamicsADTest) {
//  srand (421337);
//  ForwardDynamicsADTestTemplate(*this, 1, 1e-5);
//}

// -----------------------------------------------------------------------------

template<typename T>
void InverseDynamicsADTestTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec
) {
  Model ad_model      = obj.model;
  Model fdc_model      = obj.model;

  ADModel ad_d_model(ad_model);
  ADModel fdc_d_model(fdc_model);

  VectorNd & q     = obj.q;
  VectorNd & qdot  = obj.qdot;
  VectorNd & qddot = obj.qddot;
  VectorNd & tau   = obj.tau;
  unsigned ndirs = ad_model.qdot_size;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    q.setRandom();
    qdot.setRandom();
    qddot.setRandom();

    MatrixNd q_dirs 	= MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd qdot_dirs 	= MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd qddot_dirs = MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);

    std::vector<SpatialVector> f_ext (
          ad_model.mBodies.size(),
          SpatialVector::Zero());

    vector<vector<SpatialVector>> f_ext_dirs (
          obj.model.mBodies.size(),
          vector<SpatialVector>(
            ndirs,
            SpatialVector::Zero()));

    for (unsigned i = 0; i < f_ext.size(); i++) {
      f_ext[i].setRandom();
    }
    for (unsigned i = 0; i < f_ext.size(); i++) {
      for (unsigned idir = 0; idir < f_ext_dirs[i].size(); idir++) {
        f_ext_dirs[i][idir].setRandom();
      }
    }

    VectorNd tau_ref (tau);

    MatrixNd ad_tau  = MatrixNd::Random(ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd fd_tau  = MatrixNd::Random(ad_model.qdot_size, ad_model.qdot_size);

    AD::InverseDynamics(
          ad_model, ad_d_model,
          q, q_dirs,
          qdot, qdot_dirs,
          qddot, qddot_dirs,
          tau, ad_tau,
          &f_ext, &f_ext_dirs);

    FDC::InverseDynamics(
          fdc_model, &fdc_d_model,
          q, q_dirs,
          qdot, qdot_dirs,
          qddot, qddot_dirs,
          tau_ref, fd_tau,
          &f_ext, &f_ext_dirs);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fdc_model, fdc_d_model);

    CHECK_ARRAY_CLOSE (tau_ref.data(), tau.data(), tau_ref.rows(), array_close_prec);
    CHECK_ARRAY_CLOSE (fd_tau.data(), ad_tau.data(), fd_tau.cols()*fd_tau.rows(), array_close_prec);
  }
}

TEST_FIXTURE(CartPendulum, CartPendulumInverseDynamicsADTest){
  InverseDynamicsADTestTemplate<CartPendulum>(*this, 1, 1e-5);
}

TEST_FIXTURE(Arm2DofX, Arm2DofXInverseDynamicsADTest) {
  InverseDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm2DofZ, Arm2DofZInverseDynamicsADTest) {
  InverseDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm3DofXZYp, Arm3DofXZYpInverseDynamicsADTest) {
  InverseDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE(Arm3DofXZZp, Arm3DofXZYpInverseDynamicsADTest) {
  InverseDynamicsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE( Human36, Human36InverseDynamicsADTest) {
  srand (421337);
  InverseDynamicsADTestTemplate(*this, 1, 1e-5);
}

// -----------------------------------------------------------------------------

template<typename T>
void InverseDynamicsEDTestTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec
) {
  Model ad_model      = obj.model;
  Model ed_model      = obj.model;

  ADModel ad_d_model(ad_model);
  EDModel ed_d_model(ed_model);

  VectorNd & q     = obj.q;
  VectorNd & qdot  = obj.qdot;
  VectorNd & qddot = obj.qddot;
  VectorNd & tau   = obj.tau;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    q.setRandom();
    qdot.setRandom();
    qddot.setRandom();

    MatrixNd q_dirs   = MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd qdot_dirs  = MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd qddot_dirs = MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);

    VectorNd tau_ref (tau);

    MatrixNd ad_tau  = MatrixNd::Zero(ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd ed_tau  = MatrixNd::Zero(ad_model.qdot_size, ad_model.qdot_size);

    RigidBodyDynamics::AD::InverseDynamics(
          ad_model, ad_d_model,
          q, q_dirs,
          qdot, qdot_dirs,
          qddot, qddot_dirs,
          tau, ad_tau,
          NULL, NULL
    );

    RigidBodyDynamics::ED::InverseDynamics(
          ed_model, ed_d_model,
          q, q_dirs,
          qdot, qdot_dirs,
          qddot, qddot_dirs,
          tau_ref, ed_tau,
          NULL, NULL
    );

    for (unsigned int i = 0; i < obj.model.mBodies.size(); i++) {
      SpatialDirection v1 = ed_d_model.v[i];
      for (unsigned int idir = 0; idir < ad_model.qdot_size; idir++) {
        SpatialVector v2 = ad_d_model.v[i][idir];
        if( (v1.col(idir) - v2).norm() > 1e-6) {
          std::cout << "v " << i << "," << idir << std::endl;
          std::cout << v1.col(idir).transpose() << std::endl;
          std::cout << v2.transpose() << std::endl;
        }
        CHECK_ARRAY_CLOSE (v1.col(idir).data(), v2.data(), 6, 1e-6);
      }
    }

    for (unsigned int i = 0; i < obj.model.mBodies.size(); i++) {
      SpatialDirection c1 = ed_d_model.c[i];
      for (unsigned int idir = 0; idir < ad_model.qdot_size; idir++) {
        SpatialVector c2 = ad_d_model.c[i][idir];
        if( (c1.col(idir) - c2).norm() > 1e-6) {
          std::cout << "c " << i << "," << idir << std::endl;
          std::cout << c1.col(idir).transpose() << std::endl;
          std::cout << c2.transpose() << std::endl;
        }
        CHECK_ARRAY_CLOSE (c1.col(idir).data(), c2.data(), 6, 1e-6);
      }
    }

    for (unsigned int i = 0; i < obj.model.mBodies.size(); i++) {
      SpatialDirection a1 = ed_d_model.a[i];
      for (unsigned int idir = 0; idir < ad_model.qdot_size; idir++) {
        SpatialVector a2 = ad_d_model.a[i][idir];
        if( (a1.col(idir) - a2).norm() > 1e-6) {
          std::cout << "a " << i << "," << idir << std::endl;
          std::cout << a1.col(idir).transpose() << std::endl;
          std::cout << a2.transpose() << std::endl;
        }
        CHECK_ARRAY_CLOSE (a1.col(idir).data(), a2.data(), 6, 1e-6);
      }
    }

//    for (unsigned int i = 0; i < obj.model.mBodies.size(); i++) {
//      SpatialDirection h1 = ed_d_model.h_q[i] + ed_d_model.h_qdot[i];
//      for (unsigned int idir = 0; idir < ad_model.qdot_size; idir++) {
//        SpatialVector h2 = ad_d_model.h[i][idir];
//        if( (h1.col(idir) - h2).norm() > 1e-6) {
//          std::cout << "h " << i << "," << idir << std::endl;
//          std::cout << h1.col(idir).transpose() << std::endl;
//          std::cout << h2.transpose() << std::endl;
//        }
//        CHECK_ARRAY_CLOSE (h1.col(idir).data(), h2.data(), 6, 1e-6);
//      }
//    }

    for (unsigned int i = 0; i < obj.model.mBodies.size(); i++) {
      SpatialDirection f1 = ed_d_model.f[i];//ed_d_model.f_q[i] + ed_d_model.f_qdot[i] + ed_d_model.f_qddot[i];
      for (unsigned int idir = 0; idir < ad_model.qdot_size; idir++) {
        SpatialVector f2 = ad_d_model.f[i][idir];
        if( (f1.col(idir) - f2).norm() > 1e-6) {
          std::cout << "f " << i << "," << idir << std::endl;
          std::cout << f1.col(idir).transpose() << std::endl;
          std::cout << f2.transpose() << std::endl;
        }
        CHECK_ARRAY_CLOSE (f1.col(idir).data(), f2.data(), 6, 1e-6);
      }
    }
//     const unsigned ndirs = ad_model.qdot_size;
//     checkModelsADvsED(ndirs, ad_model, ad_d_model, ed_model, ed_d_model);

    // std::cout << "tau_nom =    " << tau.transpose() << std::endl;
    // std::cout << "tau_ref =    " << tau_ref.transpose() << std::endl;
    // std::cout << "error, ad (max):  " << (tau - tau_ref).cwiseAbs().transpose()
    //   << " (" << (tau - tau_ref).cwiseAbs().maxCoeff() << ")" << endl;
    // std::cout << endl;

    // std::cout << "ad_tau = \n" << ad_tau << std::endl;
    // std::cout << "ed_tau = \n" << ed_tau << std::endl;
    // std::cout << "error (max): \n" << (ad_tau - ed_tau).cwiseAbs()
    //   << " (" << (ad_tau - ed_tau).cwiseAbs().maxCoeff() << ")" << endl;
    // std::cout << endl;


    CHECK_ARRAY_CLOSE (tau_ref.data(), tau.data(), tau_ref.rows(), array_close_prec);
    CHECK_ARRAY_CLOSE (ed_tau.data(), ad_tau.data(), ed_tau.cols()*ed_tau.rows(), array_close_prec);
  }
}

TEST_FIXTURE(CartPendulum, CartPendulumInverseDynamicsEDTest){
  srand (421337);
  InverseDynamicsEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE(Arm2DofX, Arm2DofXInverseDynamicsEDTest) {
  srand (421337);
  InverseDynamicsEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE(Arm2DofZ, Arm2DofZInverseDynamicsEDTest) {
  srand (421337);
  InverseDynamicsEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE(Arm3DofXZYp, Arm3DofXZYpInverseDynamicsEDTest) {
  srand (421337);
  InverseDynamicsEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE(Arm3DofXZZp, Arm3DofXZZpInverseDynamicsEDTest) {
  srand (421337);
  InverseDynamicsEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE( Human36, Human36InverseDynamicsEDTest) {
  srand (421337);
  InverseDynamicsEDTestTemplate(*this, 1, 1e-5);
}

// -----------------------------------------------------------------------------

template<typename T>
void NonlinearEffectsADTestTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec
) {
  Model ad_model      = obj.model;
  Model fdc_model      = obj.model;

  ADModel ad_d_model(ad_model);
  ADModel fdc_d_model(fdc_model);

  VectorNd q(obj.model.dof_count);
  VectorNd qdot(obj.model.dof_count);
  VectorNd tau(obj.model.dof_count);

  unsigned ndirs = ad_model.qdot_size;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    q.setRandom();
    qdot.setRandom();

    MatrixNd q_dirs     = MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd qdot_dirs  = MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);

    VectorNd tau_nom (tau);
    VectorNd ad_d_tau_nom (tau);
    VectorNd fd_d_tau_nom (tau);
    MatrixNd ad_tau_der = MatrixNd::Zero(ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd fd_tau_der = MatrixNd::Zero(ad_model.qdot_size, ad_model.qdot_size);

    NonlinearEffects(ad_model, q, qdot, tau_nom);

    AD::NonlinearEffects(ad_model, ad_d_model, q, q_dirs, qdot, qdot_dirs,
                         ad_d_tau_nom, ad_tau_der);

    FDC::NonlinearEffects(fdc_model, &fdc_d_model, q, q_dirs, qdot, qdot_dirs,
                         fd_d_tau_nom, fd_tau_der);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fdc_model, fdc_d_model);

    CHECK_ARRAY_CLOSE (tau_nom.data(), ad_d_tau_nom.data(), tau_nom.rows(), array_close_prec);
    CHECK_ARRAY_CLOSE (tau_nom.data(), fd_d_tau_nom.data(), tau_nom.rows(), array_close_prec);

    CHECK_ARRAY_CLOSE (fd_tau_der.data(), ad_tau_der.data(),
                       fd_tau_der.cols() * fd_tau_der.rows(), array_close_prec);
  }
}

TEST_FIXTURE( CartPendulum, CartPendulumNonlinearEffectsADTest) {
  NonlinearEffectsADTestTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE( Arm2DofX, Arm2DofXNonlinearEffectsADTest) {
  NonlinearEffectsADTestTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE( Arm2DofZ, Arm2DofZNonlinearEffectsADTest) {
  NonlinearEffectsADTestTemplate(*this, 10, 1e-6);
}

TEST_FIXTURE( Arm3DofXZYp, Arm3DofXZYpNonlinearEffectsADTest) {
  NonlinearEffectsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE( Arm3DofXZZp, Arm3DofXZZpNonlinearEffectsADTest) {
  NonlinearEffectsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE( Human36, Human36NonlinearEffectsADTest) {
  NonlinearEffectsADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFNonlinearEffectsADTest) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  ADConstraintSet ad_constraint_set = ADConstraintSet(constraint_set, model.dof_count);
  NonlinearEffectsADTestTemplate(*this, 10, 1e-5);
}

// -----------------------------------------------------------------------------

template<typename T>
void NonlinearEffectsEDTestTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec
) {
  Model model = obj.model;
  Model ad_model = obj.model;
  Model ed_model = obj.model;

  ADModel ad_d_model(model);
  EDModel ed_d_model(model);

  VectorNd q(obj.model.dof_count);
  VectorNd qdot(obj.model.dof_count);
  VectorNd tau(obj.model.dof_count);

  // unsigned ndirs = ad_model.qdot_size;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    q.setRandom();
    qdot.setRandom();

    MatrixNd q_dirs     = MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd qdot_dirs  = MatrixNd::Random (ad_model.qdot_size, ad_model.qdot_size);

    VectorNd tau_nom (tau);
    VectorNd ad_d_tau_nom (tau);
    VectorNd ed_d_tau_nom (tau);
    MatrixNd ad_tau_der = MatrixNd::Zero(ad_model.qdot_size, ad_model.qdot_size);
    MatrixNd ed_tau_der = MatrixNd::Zero(ad_model.qdot_size, ad_model.qdot_size);

    NonlinearEffects(ad_model, q, qdot, tau_nom);

    RigidBodyDynamics::AD::NonlinearEffects(
      ad_model, ad_d_model, q, q_dirs, qdot, qdot_dirs,
      ad_d_tau_nom, ad_tau_der
    );

    RigidBodyDynamics::ED::NonlinearEffects(
      ed_model, ed_d_model, q, q_dirs, qdot, qdot_dirs,
      ed_d_tau_nom, ed_tau_der
    );

//    unsigned ndirs = ad_model.qdot_size;
//    checkModelsADvsED(ndirs, ad_model, ad_d_model, ed_model, ed_d_model);

    for (unsigned int i = 0; i < model.mBodies.size(); i++) {
      SpatialDirection v1 = ed_d_model.v[i];
      for (unsigned int idir = 0; idir < ad_model.qdot_size; idir++) {
        SpatialVector v2 = ad_d_model.v[i][idir];
        if( (v1.col(idir) - v2).norm() > 1e-6) {
          std::cout << v1.col(idir).transpose() << std::endl;
          std::cout << v2.transpose() << std::endl;
        }
        CHECK_ARRAY_CLOSE (v1.col(idir).data(), v2.data(), 6, 1e-6);
      }
    }

    for (unsigned int i = 0; i < model.mBodies.size(); i++) {
      SpatialDirection c1 = ed_d_model.c[i];
      for (unsigned int idir = 0; idir < ad_model.qdot_size; idir++) {
        SpatialVector c2 = ad_d_model.c[i][idir];
        if( (c1.col(idir) - c2).norm() > 1e-6) {
          std::cout << c1.col(idir).transpose() << std::endl;
          std::cout << c2.transpose() << std::endl;
        }
        CHECK_ARRAY_CLOSE (c1.col(idir).data(), c2.data(), 6, 1e-6);
      }
    }

    for (unsigned int i = 0; i < model.mBodies.size(); i++) {
      SpatialDirection a1 = ed_d_model.a[i];
      for (unsigned int idir = 0; idir < ad_model.qdot_size; idir++) {
        SpatialVector a2 = ad_d_model.a[i][idir];
        if( (a1.col(idir) - a2).norm() > 1e-6) {
          std::cout << a1.col(idir).transpose() << std::endl;
          std::cout << a2.transpose() << std::endl;
        }
        CHECK_ARRAY_CLOSE (a1.col(idir).data(), a2.data(), 6, 1e-6);
      }
    }

//    for (unsigned int i = 0; i < model.mBodies.size(); i++) {
//      SpatialDirection h1 = ed_d_model.h_q[i] + ed_d_model.h_qdot[i];
//      for (unsigned int idir = 0; idir < ad_model.qdot_size; idir++) {
//        SpatialVector h2 = ad_d_model.h[i][idir];
//        if( (h1.col(idir) - h2).norm() > 1e-6) {
//          std::cout << h1.col(idir).transpose() << std::endl;
//          std::cout << h2.transpose() << std::endl;
//        }
//        CHECK_ARRAY_CLOSE (h1.col(idir).data(), h2.data(), 6, 1e-6);
//      }
//    }

    for (unsigned int i = 0; i < model.mBodies.size(); i++) {
      SpatialDirection f1 = ed_d_model.f[i];
      for (unsigned int idir = 0; idir < ad_model.qdot_size; idir++) {
        SpatialVector f2 = ad_d_model.f[i][idir];
        if( (f1.col(idir) - f2).norm() > 1e-6) {
          std::cout << "f " << i << "," << idir << std::endl;
          std::cout << f1.col(idir).transpose() << std::endl;
          std::cout << f2.transpose() << std::endl;
        }
        CHECK_ARRAY_CLOSE (f1.col(idir).data(), f2.data(), 6, 1e-6);
      }
    }

    // std::cout << "tau_nom =    " << tau_nom.transpose() << std::endl;
    // std::cout << "ad_tau_nom = " << ad_d_tau_nom.transpose() << std::endl;
    // std::cout << "ed_tau_nom = " << ed_d_tau_nom.transpose() << std::endl;
    // std::cout << "error, ad (max):  " << (ad_d_tau_nom - tau_nom).cwiseAbs().transpose()
    //   << " (" << (ad_d_tau_nom - tau_nom).cwiseAbs().maxCoeff() << ")" << endl;
    // std::cout << "error, ed (max):  " << (ed_d_tau_nom - tau_nom).cwiseAbs().transpose()
    //   << " (" << (ed_d_tau_nom - tau_nom).cwiseAbs().maxCoeff() << ")" << endl;
    // std::cout << endl;

    // std::cout << "ad_tau_der = \n" << ad_tau_der << std::endl;
    // std::cout << "ed_tau_der = \n" << ed_tau_der << std::endl;
    // std::cout << "error (max): \n" << (ad_tau_der - ed_tau_der).cwiseAbs()
    //   << " (" << (ad_tau_der - ed_tau_der).cwiseAbs().maxCoeff() << ")" << endl;
    // std::cout << endl;

//    CHECK_ARRAY_CLOSE (tau_nom.data(), ad_d_tau_nom.data(), tau_nom.rows(), array_close_prec);
    CHECK_ARRAY_CLOSE (ad_d_tau_nom.data(), ed_d_tau_nom.data(), tau_nom.rows(), array_close_prec);

    CHECK_ARRAY_CLOSE (
      ed_tau_der.data(), ad_tau_der.data(),
      ed_tau_der.cols() * ed_tau_der.rows(), array_close_prec
    );
  }
}

TEST_FIXTURE( CartPendulum, CartPendulumNonlinearEffectsEDTest) {
  NonlinearEffectsEDTestTemplate(*this, 1, 1e-6);
}

 TEST_FIXTURE( Arm2DofX, Arm2DofXNonlinearEffectsEDTest) {
   NonlinearEffectsEDTestTemplate(*this, 10, 1e-6);
 }

 TEST_FIXTURE( Arm2DofZ, Arm2DofZNonlinearEffectsEDTest) {
   NonlinearEffectsEDTestTemplate(*this, 10, 1e-6);
 }

 TEST_FIXTURE( Arm3DofXZYp, Arm3DofXZYpNonlinearEffectsEDTest) {
   NonlinearEffectsEDTestTemplate(*this, 10, 1e-5);
 }

 TEST_FIXTURE( Arm3DofXZZp, Arm3DofXZZpNonlinearEffectsEDTest) {
   NonlinearEffectsEDTestTemplate(*this, 10, 1e-5);
 }

 TEST_FIXTURE( Human36, Human36NonlinearEffectsEDTest) {
   NonlinearEffectsEDTestTemplate(*this, 10, 1e-5);
 }

 TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFNonlinearEffectsEDTest) {
   // add contacts and bind them to constraint set
   constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
   constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
   constraint_set.Bind (model);
   EDConstraintSet ad_constraint_set(constraint_set, model.dof_count);
   NonlinearEffectsEDTestTemplate(*this, 10, 1e-5);
 }

// -----------------------------------------------------------------------------
template<typename T>
void CompositeRigidBodyAlgorithmADvsFDTestTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec) {

  Model model = obj.model;
  Model fd_model = obj.model;
  Model ad_model = obj.model;
  ADModel ad_d_model(ad_model);
  ADModel fd_d_model(fd_model);
  for (unsigned trial = 0; trial < trial_count; trial++) {
    unsigned ndirs = model.qdot_size;
    // initialize inputs and directions
    VectorNd q = VectorNd::Random(model.qdot_size);
    MatrixNd q_dirs = MatrixNd::Random (model.qdot_size, ndirs);

    // initialize outputs and derivatives
    MatrixNd H = MatrixNd::Zero(q.size(),q.size());
    MatrixNd ad_H = MatrixNd::Zero(model.q_size, ndirs);
    MatrixNd fd_H = MatrixNd::Zero(model.q_size, ndirs);
    vector<MatrixNd> fd_H_dirs(ndirs, MatrixNd::Zero(model.q_size,model.q_size));
    vector<MatrixNd> ad_H_dirs(ndirs, MatrixNd::Zero(model.q_size,model.q_size));

    // compute nominal value
    CompositeRigidBodyAlgorithm(model, q, H, true);

    // evaluate derivatives
    // -> using finite differences
    FDC::CompositeRigidBodyAlgorithm(fd_model, &fd_d_model, q, q_dirs, fd_H, fd_H_dirs);

    // -> using algorithmic differences
    AD::CompositeRigidBodyAlgorithm(ad_model, ad_d_model, q, q_dirs, ad_H, ad_H_dirs);

    // check derivatives of model quantities
    checkModelsADvsFD(ndirs,
                      ad_model, ad_d_model,
                      fd_model, fd_d_model);

    // check nominal values for consistency
    CHECK_ARRAY_CLOSE (H.data(), ad_H.data(),
                       model.q_size * model.q_size,
                       1e-16);
    CHECK_ARRAY_CLOSE (H.data(), fd_H.data(),
                       model.q_size * model.q_size,
                       1e-16);

    // check AD vs FD derivatives for consistency
    for (size_t idir = 0; idir < model.qdot_size; idir++) {
      CHECK_ARRAY_CLOSE (fd_H_dirs[idir].data(), ad_H_dirs[idir].data(),
                         model.q_size * model.q_size,
                         array_close_prec);
    }
  }
}

TEST_FIXTURE( CartPendulum, CartPendulumCompositeRigidBodyAlgorithmADvsFDTest) {
  CompositeRigidBodyAlgorithmADvsFDTestTemplate(*this, 1, 1e-8);
}

TEST_FIXTURE( Arm2DofX, Arm2DofXCompositeRigidBodyAlgorithmADvsFDTest) {
  CompositeRigidBodyAlgorithmADvsFDTestTemplate(*this, 1, 1e-8);
}

TEST_FIXTURE( Arm2DofZ, Arm2DofZCompositeRigidBodyAlgorithmADvsFDTest) {
  CompositeRigidBodyAlgorithmADvsFDTestTemplate(*this, 1, 1e-8);
}

TEST_FIXTURE( Arm3DofXZYp, Arm3DofXZYpCompositeRigidBodyAlgorithmADvsFDTest) {
  CompositeRigidBodyAlgorithmADvsFDTestTemplate(*this, 1, 1e-6);
}

TEST_FIXTURE( Arm3DofXZZp, Arm3DofXZZpCompositeRigidBodyAlgorithmADvsFDTest) {
  CompositeRigidBodyAlgorithmADvsFDTestTemplate(*this, 1, 1e-6);
}

TEST_FIXTURE( FixedBase6DoF9DoF, FixedBase6DoF9DoFCompositeRigidBodyAlgorithmADvsFDTest) {
  CompositeRigidBodyAlgorithmADvsFDTestTemplate(*this, 1, 1e-6);
}

TEST_FIXTURE( Human36, Human36CompositeRigidBodyAlgorithmADTest) {
  CompositeRigidBodyAlgorithmADvsFDTestTemplate(*this, 1, 1e-5);
}


// -----------------------------------------------------------------------------
template<typename T>
void CompositeRigidBodyAlgorithmEDTestTemplate(
  T & obj,
  unsigned trial_count,
  double array_close_prec
) {
  Model model = obj.model;
  Model ed_model = obj.model;
  Model ad_model = obj.model;
  ADModel ad_d_model(model);
  EDModel ed_d_model(model);
  for (unsigned trial = 0; trial < trial_count; trial++) {
    unsigned ndirs = model.qdot_size;
    // initialize inputs and directions
    VectorNd q = VectorNd::Random(model.qdot_size);
    MatrixNd q_dirs = MatrixNd::Random (model.qdot_size, ndirs);

    // initialize outputs and derivatives
    MatrixNd H = MatrixNd::Zero(q.size(),q.size());
    MatrixNd ad_H = MatrixNd::Zero(model.q_size, ndirs);
    MatrixNd ed_H = MatrixNd::Zero(model.q_size, ndirs);
    vector<MatrixNd> ed_H_dirs(ndirs, MatrixNd::Zero(model.q_size,model.q_size));
    vector<MatrixNd> ad_H_dirs(ndirs, MatrixNd::Zero(model.q_size,model.q_size));

    // compute nominal value
    CompositeRigidBodyAlgorithm(model, q, H, true);

    // evaluate derivatives
    // -> using finite differences
    RigidBodyDynamics::ED::CompositeRigidBodyAlgorithm(
      ed_model, ed_d_model, q, q_dirs, ed_H, ed_H_dirs
    );

    // -> using algorithmic differences
    RigidBodyDynamics::AD::CompositeRigidBodyAlgorithm(
      ad_model, ad_d_model, q, q_dirs, ad_H, ad_H_dirs
    );

    //     check derivatives of model quantities
//     checkModelsADvsED(
//       ndirs,
//       ad_model, ad_d_model,
//       ed_model, ed_d_model
//     );

    // check nominal values for consistency
    const double NOM_TOL = 1e-16;
    MatrixNd error;
    double max;

    CHECK_ARRAY_CLOSE (
      H.data(), ad_H.data(),
                       model.q_size * model.q_size,
      NOM_TOL
    );

    error = (H - ad_H).cwiseAbs();
    max = error.maxCoeff();
    if (max > NOM_TOL) {
      std::cout << "error = \n" << error << std::endl;
      std::cout << "max   = " << max << std::endl;
      std::cout << "H = \n"<< H << std::endl;
      std::cout << "ad_H = \n"<< ad_H << std::endl;
      std::cout << endl;
    }
    error = (H - ed_H).cwiseAbs();
    max = error.maxCoeff();
    if (max > NOM_TOL) {
      std::cout << "error = \n" << error << std::endl;
      std::cout << "max   = " << max << std::endl;
      std::cout << "H = \n"<< H << std::endl;
      std::cout << "ed_H = \n"<< ed_H << std::endl;
      std::cout << endl;
    }
    CHECK_ARRAY_CLOSE (
      H.data(), ed_H.data(),
      model.q_size * model.q_size,
      NOM_TOL
    );

    error = (ad_H - ed_H).cwiseAbs();
    max = error.maxCoeff();
    if (max > NOM_TOL) {
      std::cout << "error = \n" << error << std::endl;
      std::cout << "max   = " << max << std::endl;
      std::cout << "ad_H = \n"<< ad_H << std::endl;
      std::cout << "ed_H = \n"<< ed_H << std::endl;
      std::cout << endl;
    }

    CHECK_ARRAY_CLOSE (
      ad_H.data(), ed_H.data(),
      model.q_size * model.q_size,
      NOM_TOL
    );

    // check AD vs ED derivatives for consistency
    for (size_t idir = 0; idir < model.qdot_size; idir++) {
      error = (ed_H_dirs[idir] - ad_H_dirs[idir]).cwiseAbs();
      max = error.maxCoeff();
      if (max > 1e-12) {
        std::cout << "error [" << idir << "]= \n" << error << std::endl;
        std::cout << "max   [" << idir << "]= " << max << std::endl;
        // std::cout << "ed_H_dirs[" << idir << "]=\n"
        //   << ed_H_dirs[idir] << std::endl;
        // std::cout << "ad_H_dirs[" << idir << "]=\n"
        //   << ad_H_dirs[idir] << std::endl;
        std::cout << endl;
      }
      CHECK_ARRAY_CLOSE (
        ed_H_dirs[idir].data(), ad_H_dirs[idir].data(),
                         model.q_size * model.q_size,
        array_close_prec
      );
    }
  }
}


TEST_FIXTURE (CartPendulum, CartPendulumCompositeRigidBodyAlgorithmEDTest)
{
  CompositeRigidBodyAlgorithmEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE (Arm2DofX, Arm2DofXCompositeRigidBodyAlgorithmEDTest)
{
  CompositeRigidBodyAlgorithmEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE (Arm2DofZ, Arm2DofZCompositeRigidBodyAlgorithmEDTest)
{
  CompositeRigidBodyAlgorithmEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE (Arm3DofXZYp, Arm3DofXZYpCompositeRigidBodyAlgorithmEDTest)
{
  CompositeRigidBodyAlgorithmEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE (Arm3DofXZZp, Arm3DofXZZpCompositeRigidBodyAlgorithmEDTest)
{
  CompositeRigidBodyAlgorithmEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE (
  FixedBase6DoF9DoF,
  FixedBase6DoF9DoFCompositeRigidBodyAlgorithmEDTest
) {
  CompositeRigidBodyAlgorithmEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE(
  MultiPendulumWithBranches,
  MultiPendulumWithBranchesCompositeRigidBodyAlgorithmEDTest_4_2
) {
  this->create_model(4, 2);
  CompositeRigidBodyAlgorithmEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE(
  MultiPendulumWithBranches,
  MultiPendulumWithBranchesCompositeRigidBodyAlgorithmEDTest_6_3
) {
  this->create_model(6, 3);
  CompositeRigidBodyAlgorithmEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE(
  MultiPendulumWithBranches,
  MultiPendulumWithBranchesCompositeRigidBodyAlgorithmEDTest_12_3
) {
  this->create_model(12, 3);
  CompositeRigidBodyAlgorithmEDTestTemplate(*this, 1, 1e-10);
}

TEST_FIXTURE( Human36, Human36CompositeRigidBodyAlgorithmEDTest) {
  CompositeRigidBodyAlgorithmEDTestTemplate(*this, 1, 1e-10);
}

// -----------------------------------------------------------------------------
