#include <UnitTest++.h>

#include <iostream>
#include <iomanip>

#include "rbdl/Model.h"
#include "rbdl/Kinematics.h"
#include "rbdl/rbdl_mathutils.h"

#include "ConstraintsED.h"
#include "ConstraintsAD.h"
#include "ConstraintsFD.h"
#include "ConstraintsFDC.h"


#include "DynamicsAD.h"

#include "ModelCheckADvsFD.h"

#include "Fixtures.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

double const TEST_PREC = 1.0e-8;

// -----------------------------------------------------------------------------

template <typename T>
void CalcContactJacobianTemplate(
    T & obj,
    const unsigned trial_count,
    const double array_close_prec
) {
  Model   ad_model   = obj.model;
  Model   fd_model   = obj.model;

  ADModel ad_d_model = ADModel(ad_model);
  ADModel fd_d_model = ADModel(fd_model);

  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet fd_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = ADConstraintSet(ad_cs, ad_model.dof_count);
  ADConstraintSet fd_d_cs = ADConstraintSet(fd_cs, fd_model.dof_count);

  // set up input quantities
  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;
  bool update_kinematics = true;

  VectorNd q = VectorNd::Zero(nq);
  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  for (unsigned trial = 0; trial < trial_count; trial++) {
    // set up no output quantities
    MatrixNd G = MatrixNd::Zero (3, nq);
    MatrixNd G_ad = MatrixNd::Zero (3, nq);
    MatrixNd G_fd = MatrixNd::Zero (3, nq);

    // set up derivative output quantities
    vector<MatrixNd> ad_d_G (ndirs, G_ad);
    vector<MatrixNd> fd_d_G (ndirs, G_fd);

    // call nominal version
    CalcConstraintsJacobian(ad_model, q, ad_cs, G, update_kinematics);

    AD::CalcConstraintsJacobian(
          ad_model, ad_d_model,
          q, q_dirs,
          ad_cs, ad_d_cs,
          G_ad, ad_d_G,
          update_kinematics
          );

    FDC::CalcConstraintsJacobian(
          fd_model, &fd_d_model,
          q, q_dirs,
          fd_cs, &fd_d_cs,
          G_fd, fd_d_G
          );

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);
    checkConstraintSetsADvsFD(ndirs, ad_cs, ad_d_cs, fd_cs, fd_d_cs);

    CHECK_ARRAY_CLOSE(G.data(), G_ad.data(), G.size(), array_close_prec);
    CHECK_ARRAY_CLOSE(G.data(), G_fd.data(), G.size(), array_close_prec);
    CHECK_ARRAY_CLOSE(G_fd.data(), G_ad.data(), G.size(), array_close_prec);

    for (unsigned idir = 0; idir < ndirs; idir++) {
      CHECK_ARRAY_CLOSE(
            fd_d_G[idir].data(),
            ad_d_G[idir].data(),
            ad_d_G[idir].size(),
            array_close_prec);
    }
    q.setRandom();
    q_dirs.setRandom();
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFCalcContactJacobian) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  CalcContactJacobianTemplate(*this, 10, 1e-9);
}

// -----------------------------------------------------------------------------


template <typename T>
void CalcContactJacobianEDvsADTemplate(
    T & obj,
    const unsigned trial_count,
    const double array_close_prec
) {
  Model   ad_model   = obj.model;
  Model   ed_model   = obj.model;

  ADModel ad_d_model = ADModel(ad_model);
  EDModel ed_d_model = EDModel(ed_model);

  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet ed_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = ADConstraintSet(ad_cs, ad_model.dof_count);
  EDConstraintSet ed_d_cs = EDConstraintSet(ed_cs, ed_model.dof_count);

  // set up input quantities
  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;
  bool update_kinematics = true;

  VectorNd q = VectorNd::Zero(nq);
  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  for (unsigned trial = 0; trial < trial_count; trial++) {
    // set up no output quantities
    MatrixNd G = MatrixNd::Zero (3, nq);
    MatrixNd G_ad = MatrixNd::Zero (3, nq);
    MatrixNd G_fd = MatrixNd::Zero (3, nq);

    // set up derivative output quantities
    vector<MatrixNd> ad_d_G (ndirs, G_ad);
    vector<MatrixNd> ed_d_G (ndirs, G_fd);

    // call nominal version
    CalcConstraintsJacobian(ad_model, q, ad_cs, G, update_kinematics);

    RigidBodyDynamics::AD::CalcConstraintsJacobian(
          ad_model, ad_d_model,
          q, q_dirs,
          ad_cs, ad_d_cs,
          G_ad, ad_d_G,
          update_kinematics
          );

    RigidBodyDynamics::ED::CalcConstraintsJacobian(
          ed_model, ed_d_model,
          q, q_dirs,
          ed_cs, ed_d_cs,
          G_fd, ed_d_G
          );

    checkModelsADvsED(ndirs, ad_model, ad_d_model, ed_model, ed_d_model);
    checkConstraintSetsADvsED(ndirs, ad_cs, ad_d_cs, ed_cs, ed_d_cs);

    CHECK_ARRAY_CLOSE(G.data(), G_ad.data(), G.size(), array_close_prec);
    CHECK_ARRAY_CLOSE(G.data(), G_fd.data(), G.size(), array_close_prec);
    CHECK_ARRAY_CLOSE(G_fd.data(), G_ad.data(), G.size(), array_close_prec);

    for (unsigned idir = 0; idir < ndirs; idir++) {
      CHECK_ARRAY_CLOSE(
            ed_d_G[idir].data(),
            ad_d_G[idir].data(),
            ad_d_G[idir].size(),
            array_close_prec);
    }
    q.setRandom();
    q_dirs.setRandom();
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFCalcContactJacobianEDvsAD) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  CalcContactJacobianEDvsADTemplate(*this, 10, 1e-12);
}

// -----------------------------------------------------------------------------


template <typename T>
void CalcContactSystemVariablesTemplate(
    T & obj,
    unsigned trial_count) {
  Model   model      = obj.model;
  Model   ad_model   = obj.model;
  Model   fd_model   = obj.model;

  ADModel ad_d_model = ADModel(ad_model);
  ADModel fd_d_model = ADModel(fd_model);

  ConstraintSet cs = obj.constraint_set;
  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet fd_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = ADConstraintSet(ad_cs, ad_model.dof_count);
  ADConstraintSet fd_d_cs = ADConstraintSet(fd_cs, fd_model.dof_count);

  // set up input quantities
  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;

  VectorNd q = VectorNd::Zero(nq);
  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  VectorNd qd = VectorNd::Zero(nq);
  MatrixNd qd_dirs = MatrixNd::Identity(nq, nq);
  VectorNd tau = VectorNd::Zero(nq);
  MatrixNd tau_dirs = MatrixNd::Identity(nq, nq);
  for (unsigned trial = 0; trial < trial_count; trial++) {
    CalcConstrainedSystemVariables(model, q, qd, tau, cs);

    AD::CalcConstrainedSystemVariables(
      ad_model, ad_d_model,
      q, q_dirs, qd, qd_dirs,
      tau, tau_dirs, ad_cs, ad_d_cs
    );

    FDC::CalcConstrainedSystemVariables(
      fd_model, &fd_d_model,
      q, q_dirs, qd, qd_dirs,
      tau, tau_dirs, fd_cs, &fd_d_cs
    );

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);
    checkConstraintSetsADvsFD(ndirs, ad_cs, ad_d_cs, fd_cs, fd_d_cs);

    q.setRandom();
    qd.setRandom();
    tau.setRandom();
    q_dirs.setRandom();
    qd_dirs.setRandom();
    tau_dirs.setRandom();
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFCalcContactSystemVariables) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  CalcContactSystemVariablesTemplate(*this, 10);
}


// -----------------------------------------------------------------------------

template <typename T>
void CalcContactSystemVariablesEDvsADTemplate(
    T & obj,
    const unsigned trial_count
) {
  Model model = obj.model;
  Model ad_model = obj.model;
  Model fd_model = obj.model;

  ADModel ad_d_model = ADModel(ad_model);
  EDModel fd_d_model = EDModel(fd_model);

  ConstraintSet cs = obj.constraint_set;
  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet fd_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = ADConstraintSet(ad_cs, ad_model.dof_count);
  EDConstraintSet fd_d_cs = EDConstraintSet(fd_cs, fd_model.dof_count);

  // set up input quantities
  int const nq = ad_model.dof_count;
  // unsigned const ndirs = nq;

  VectorNd q = VectorNd::Zero(nq);
  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  VectorNd qd = VectorNd::Zero(nq);
  MatrixNd qd_dirs = MatrixNd::Identity(nq, nq);
  VectorNd tau = VectorNd::Zero(nq);
  MatrixNd tau_dirs = MatrixNd::Identity(nq, nq);

  for (unsigned trial = 0; trial < trial_count; trial++)
  {
    CalcConstrainedSystemVariables(model, q, qd, tau, cs);

    RigidBodyDynamics::AD::CalcConstrainedSystemVariables (
      ad_model, ad_d_model,
      q, q_dirs, qd, qd_dirs,
      tau, tau_dirs, ad_cs, ad_d_cs
    );

    RigidBodyDynamics::ED::CalcConstrainedSystemVariables(
      fd_model, fd_d_model,
      q, q_dirs, qd, qd_dirs,
      tau, tau_dirs, fd_cs, fd_d_cs
    );

    // checkModelsADvsED(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);
    // checkConstraintSetsADvsED(ndirs, ad_cs, ad_d_cs, fd_cs, fd_d_cs);

    q.setRandom();
    qd.setRandom();
    tau.setRandom();
    q_dirs.setRandom();
    qd_dirs.setRandom();
    tau_dirs.setRandom();
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFCalcContactSystemVariablesEDvsAD) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  CalcContactSystemVariablesEDvsADTemplate(*this, 10);
}


// -----------------------------------------------------------------------------

template <typename T>
void ForwardDynamicsConstraintsDirectEDvsADTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec
) {
  Model   model    = obj.model;
  Model   ad_model = obj.model;
  Model   ed_model = obj.model;

  ADModel ad_d_model = ADModel(ad_model);
  EDModel ed_d_model = EDModel(ed_model);

  ConstraintSet cs = obj.constraint_set;
  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet ed_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = ADConstraintSet(ad_cs, ad_model.dof_count);
  EDConstraintSet ed_d_cs = EDConstraintSet(ed_cs, ed_model.dof_count);

  // set up input quantities
  int const nq = ad_model.dof_count;
  // unsigned const ndirs = nq;

  VectorNd qdd    = VectorNd::Zero(nq);
  VectorNd ad_qdd = VectorNd::Zero(nq);
  VectorNd ed_qdd = VectorNd::Zero(nq);
  MatrixNd ad_qdd_dirs = MatrixNd::Zero(nq, nq);
  MatrixNd ed_qdd_dirs = MatrixNd::Zero(nq, nq);

  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd qd_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd tau_dirs = MatrixNd::Identity(nq, nq);

  for (unsigned trial = 0; trial < trial_count; trial++) {
    VectorNd q = VectorNd::Random(nq);
    VectorNd qd = VectorNd::Random(nq);
    VectorNd tau = VectorNd::Random(nq);

    RigidBodyDynamics::ForwardDynamicsConstraintsDirect(model, q, qd, tau, cs, qdd);

    RigidBodyDynamics::AD::ForwardDynamicsConstraintsDirect(
          ad_model, ad_d_model,
          q, q_dirs,
          qd, qd_dirs,
          tau, tau_dirs,
          ad_cs, ad_d_cs,
          ad_qdd, ad_qdd_dirs);

    RigidBodyDynamics::ED::ForwardDynamicsConstraintsDirect(
          ed_model, ed_d_model,
          q, q_dirs,
          qd, qd_dirs,
          tau, tau_dirs,
          ed_cs, ed_d_cs,
          ed_qdd, ed_qdd_dirs);

    // checkModelsADvsED(ndirs, ad_model, ad_d_model, ed_model, ed_d_model);
    // checkConstraintSetsADvsED(ndirs, ad_cs, ad_d_cs, ed_cs, ed_d_cs);

    CHECK_ARRAY_CLOSE(qdd.data(), ad_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(qdd.data(), ed_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(ad_qdd_dirs.data(), ed_qdd_dirs.data(), nq*nq,
                      array_close_prec);

    q_dirs.setRandom();
    qd_dirs.setRandom();
    tau_dirs.setRandom();
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFForwardDynamicsConstraintsDirectEDvsAD) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  ForwardDynamicsConstraintsDirectEDvsADTemplate(*this, 1, 1e-12);
}
// -----------------------------------------------------------------------------

template <typename T>
void ForwardDynamicsConstraintsDirectTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec
) {
  Model   model    = obj.model;
  Model   ad_model = obj.model;
  Model   fd_model = obj.model;

  ADModel ad_d_model = ADModel(ad_model);
  ADModel fd_d_model = ADModel(fd_model);

  ConstraintSet cs = obj.constraint_set;
  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet fd_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = ADConstraintSet(ad_cs, ad_model.dof_count);
  ADConstraintSet fd_d_cs = ADConstraintSet(fd_cs, fd_model.dof_count);

  // set up input quantities
  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;

  VectorNd qdd    = VectorNd::Zero(nq);
  VectorNd ad_qdd = VectorNd::Zero(nq);
  VectorNd fd_qdd = VectorNd::Zero(nq);
  MatrixNd ad_qdd_dirs = MatrixNd::Zero(nq, nq);
  MatrixNd fd_qdd_dirs = MatrixNd::Zero(nq, nq);

  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd qd_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd tau_dirs = MatrixNd::Identity(nq, nq);

  for (unsigned trial = 0; trial < trial_count; trial++) {
    VectorNd q = VectorNd::Random(nq);
    VectorNd qd = VectorNd::Random(nq);
    VectorNd tau = VectorNd::Random(nq);

    RigidBodyDynamics::ForwardDynamicsConstraintsDirect(model, q, qd, tau, cs, qdd);

    RigidBodyDynamics::AD::ForwardDynamicsConstraintsDirect(
          ad_model, ad_d_model,
          q, q_dirs,
          qd, qd_dirs,
          tau, tau_dirs,
          ad_cs, ad_d_cs,
          ad_qdd, ad_qdd_dirs);

    RigidBodyDynamics::FDC::ForwardDynamicsConstraintsDirect(
          fd_model, &fd_d_model,
          q, q_dirs,
          qd, qd_dirs,
          tau, tau_dirs,
          fd_cs, &fd_d_cs,
          fd_qdd, fd_qdd_dirs);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);
    checkConstraintSetsADvsFD(ndirs, ad_cs, ad_d_cs, fd_cs, fd_d_cs);

    CHECK_ARRAY_CLOSE(qdd.data(), ad_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(qdd.data(), fd_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(ad_qdd_dirs.data(), fd_qdd_dirs.data(), nq*nq,
                      array_close_prec);

    q_dirs.setRandom();
    qd_dirs.setRandom();
    tau_dirs.setRandom();
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFForwardDynamicsConstraintsDirect) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  ForwardDynamicsConstraintsDirectTemplate(*this, 10, 1e-7);
}

// -----------------------------------------------------------------------------

template <typename T>
void SolveConstrainedSystemDirectTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec
) {
  Model   model    = obj.model;
  Model   ad_model = obj.model;
  Model   fd_model = obj.model;

  ADModel ad_d_model = ADModel(ad_model);
  ADModel fd_d_model = ADModel(fd_model);

  ConstraintSet cs = obj.constraint_set;
  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet fd_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = ADConstraintSet(ad_cs, ad_model.dof_count);
  ADConstraintSet fd_d_cs = ADConstraintSet(fd_cs, fd_model.dof_count);

  // set up input quantities
  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;

  VectorNd qdd    = VectorNd::Zero(nq);
  VectorNd ad_qdd = VectorNd::Zero(nq);
  VectorNd fd_qdd = VectorNd::Zero(nq);
  MatrixNd ad_qdd_dirs = MatrixNd::Zero(nq, nq);
  MatrixNd fd_qdd_dirs = MatrixNd::Zero(nq, nq);

  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd qd_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd tau_dirs = MatrixNd::Identity(nq, nq);

  for (unsigned trial = 0; trial < trial_count; trial++) {
    VectorNd q = VectorNd::Random(nq);
    VectorNd qd = VectorNd::Random(nq);
    VectorNd tau = VectorNd::Random(nq);

    // NOTE we need to compute the system variables first
    RigidBodyDynamics::CalcConstrainedSystemVariables(model, q, qd, tau, cs);

    RigidBodyDynamics::AD::CalcConstrainedSystemVariables(
      fd_model, fd_d_model,
      q, q_dirs,
      qd, qd_dirs,
      tau, tau_dirs,
      fd_cs, fd_d_cs
    );

    RigidBodyDynamics::AD::CalcConstrainedSystemVariables(
      ad_model, ad_d_model,
      q, q_dirs,
      qd, qd_dirs,
      tau, tau_dirs,
      ad_cs, ad_d_cs
    );

    VectorNd c = tau - cs.C;
    MatrixNd cdot(c.rows(), ndirs);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      cdot.col(idir) = tau_dirs.col(idir) - ad_cs.C.col(idir);
    }

    CHECK_ARRAY_CLOSE(cs.H.data(), fd_cs.H.data(), cs.H.size(), array_close_prec);
    CHECK_ARRAY_CLOSE(cs.H.data(), ad_cs.H.data(), cs.H.size(), array_close_prec);
    CHECK_ARRAY_CLOSE(cs.G.data(), fd_cs.G.data(), cs.G.size(), array_close_prec);
    CHECK_ARRAY_CLOSE(cs.G.data(), ad_cs.G.data(), cs.G.size(), array_close_prec);
    CHECK_ARRAY_CLOSE(cs.gamma.data(), fd_cs.gamma.data(), cs.gamma.size(), array_close_prec);
    CHECK_ARRAY_CLOSE(cs.gamma.data(), ad_cs.gamma.data(), cs.gamma.size(), array_close_prec);

    RigidBodyDynamics::SolveConstrainedSystemDirect(
        cs.H,
        cs.G,
        c,
        cs.gamma,
        qdd,
        cs.force,
        cs.A,
        cs.b,
        cs.x,
        cs.linear_solver
    );

    RigidBodyDynamics::AD::SolveConstrainedSystemDirect(
         ad_cs.H, ad_d_cs.H,
         ad_cs.G, ad_d_cs.G,
         c, cdot,
         ad_cs.gamma, ad_d_cs.gamma,
         ad_cs.A, ad_d_cs.A,
         ad_cs.b, ad_d_cs.b,
         ad_cs.x, ad_d_cs.x,
         ad_cs.linear_solver,
         ndirs
    );

    ad_qdd = ad_cs.x.segment(0, ad_model.dof_count);
    ad_qdd_dirs = ad_d_cs.x.block(0, 0, ad_model.dof_count, ndirs);

    RigidBodyDynamics::FDC::SolveConstrainedSystemDirect(
         fd_cs.H, fd_d_cs.H,
         fd_cs.G, fd_d_cs.G,
         c, cdot,
         fd_cs.gamma, fd_d_cs.gamma,
         fd_cs.A, fd_d_cs.A,
         fd_cs.b, fd_d_cs.b,
         fd_cs.x, fd_d_cs.x,
         fd_cs.linear_solver,
         ndirs
    );

    qdd = cs.x.segment(0, model.dof_count);
    fd_qdd = fd_cs.x.segment(0, fd_model.dof_count);
    fd_qdd_dirs = fd_d_cs.x.block(0, 0, fd_model.dof_count, ndirs);

    CHECK_ARRAY_CLOSE(qdd.data(), ad_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(ad_qdd.data(), fd_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(ad_qdd_dirs.data(), fd_qdd_dirs.data(), nq*nq,
                      array_close_prec);

    q_dirs.setRandom();
    qd_dirs.setRandom();
    tau_dirs.setRandom();
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFSolveConstrainedSystemDirect) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  SolveConstrainedSystemDirectTemplate(*this, 10, 1e-7);
}

// -----------------------------------------------------------------------------

template <typename T>
void ForwardDynamicsAccelerationDeltasTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec
) {
  Model   model    = obj.model;
  Model   ad_model = obj.model;
  Model   fd_model = obj.model;

  ADModel ad_d_model = ADModel(ad_model);
  ADModel fd_d_model = ADModel(fd_model);

  ConstraintSet cs = obj.constraint_set;
  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet fd_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = ADConstraintSet(ad_cs, ad_model.dof_count);
  ADConstraintSet fd_d_cs = ADConstraintSet(fd_cs, fd_model.dof_count);

  // set up input quantities
  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;

  VectorNd qdd    = VectorNd::Zero(nq);
  VectorNd ad_qdd = VectorNd::Zero(nq);
  VectorNd fd_qdd = VectorNd::Zero(nq);
  MatrixNd ad_qdd_dirs = MatrixNd::Zero(nq, nq);
  MatrixNd fd_qdd_dirs = MatrixNd::Zero(nq, nq);

  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd qd_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd tau_dirs = MatrixNd::Identity(nq, nq);

  for (unsigned trial = 0; trial < trial_count; trial++) {
    VectorNd q = VectorNd::Random(nq);
    VectorNd qd = VectorNd::Random(nq);
    VectorNd tau = VectorNd::Random(nq);

    q_dirs.setRandom();
    qd_dirs.setRandom();
    tau_dirs.setRandom();
    std::vector<SpatialVector> f_t
      = std::vector<SpatialVector> (model.mBodies.size(), SpatialVector::Random());
    std::vector<MatrixNd> ad_f_t
      = std::vector<MatrixNd> (model.mBodies.size(), SpatialMatrix::Random(6, ndirs));

    UpdateKinematics(model, q, qd, qdd);
    UpdateKinematics(ad_model, q, qd, qdd);
    UpdateKinematics(fd_model, q, qd, qdd);

    ForwardDynamics(model, q, qd, tau, cs.QDDot_0);
    ForwardDynamics(fd_model, q, qd, tau, cs.QDDot_0);
    ForwardDynamics(ad_model, q, qd, tau, cs.QDDot_0);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);

    ForwardDynamicsAccelerationDeltas(
          model, cs, qdd, obj.contact_body_id, f_t);

    AD::ForwardDynamicsAccelerationDeltas (
          ad_model, ad_d_model,
          ad_cs, ad_d_cs,
          ad_qdd, ad_qdd_dirs,
          obj.contact_body_id,
          f_t, ad_f_t
    );

    FD::ForwardDynamicsAccelerationDeltas (
          fd_model, &fd_d_model,
          fd_cs, fd_d_cs,
          fd_qdd, fd_qdd_dirs,
          obj.contact_body_id,
          f_t, ad_f_t
    );

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);
    checkModelsADvsFD(ndirs, ad_model, ad_d_model, model, ad_d_model);

    checkConstraintSetsADvsFD(ndirs, ad_cs, ad_d_cs, fd_cs, fd_d_cs);

    CHECK_ARRAY_CLOSE(ad_qdd.data(), qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(ad_qdd.data(), fd_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(ad_qdd_dirs.data(), fd_qdd_dirs.data(), nq*nq,
                      array_close_prec);
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFForwardDynamicsAccelerationDeltas) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  ForwardDynamicsAccelerationDeltasTemplate(*this, 1, 1e-6);
}

// -----------------------------------------------------------------------------

template <typename T>
void ForwardDynamicsApplyConstraintForcesTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec
) {
  Model   model    = obj.model;
  Model   ad_model = obj.model;
  Model   fd_model = obj.model;

  ADModel ad_d_model = ADModel(ad_model);
  ADModel fd_d_model = ADModel(fd_model);

  ConstraintSet cs = obj.constraint_set;
  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet fd_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = ADConstraintSet(ad_cs, ad_model.dof_count);
  ADConstraintSet fd_d_cs = ADConstraintSet(fd_cs, fd_model.dof_count);

  // set up input quantities
  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;

  VectorNd qdd    = VectorNd::Zero(nq);
  VectorNd ad_qdd = VectorNd::Zero(nq);
  VectorNd fd_qdd = VectorNd::Zero(nq);
  MatrixNd ad_qdd_dirs = MatrixNd::Zero(nq, nq);
  MatrixNd fd_qdd_dirs = MatrixNd::Zero(nq, nq);

  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd qd_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd tau_dirs = MatrixNd::Identity(nq, nq);

  for (unsigned trial = 0; trial < trial_count; trial++) {
    VectorNd q = VectorNd::Random(nq);
    VectorNd qd = VectorNd::Random(nq);
    VectorNd tau = VectorNd::Random(nq);
    UpdateKinematics(model, q, qd, qdd);
    UpdateKinematics(ad_model, q, qd, qdd);
    UpdateKinematics(fd_model, q, qd, qdd);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);

    ForwardDynamics(model, q, qd, tau, cs.QDDot_0);
    ForwardDynamics(fd_model, q, qd, tau, cs.QDDot_0);
    ForwardDynamics(ad_model, q, qd, tau, cs.QDDot_0);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);

    for (unsigned ci = 0; ci < cs.size(); ci++) {
      unsigned int movable_body_id = 0;
      switch (cs.constraintType[ci]) {
        case ConstraintSet::ContactConstraint:
          movable_body_id = GetMovableBodyId(model, cs.body[ci]);
          cs.f_ext_constraints[movable_body_id] = SpatialVector::Random();
          ad_cs.f_ext_constraints[movable_body_id] = cs.f_ext_constraints[movable_body_id];
          fd_cs.f_ext_constraints[movable_body_id] = cs.f_ext_constraints[movable_body_id];
        break;

        default:
          std::cerr << "Forward Dynamic Contact Kokkevis: unsupported constraint \
            type." << std::endl;
          assert(false);
          abort();
        break;
      }
    }

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);

    ForwardDynamicsApplyConstraintForces(model, tau, cs, qdd);

    AD::ForwardDynamicsApplyConstraintForces(
      ad_model, ad_d_model,
      tau, tau_dirs,
      ad_cs, ad_d_cs,
      ad_qdd, ad_qdd_dirs
    );

    FD::ForwardDynamicsApplyConstraintForces(
      fd_model, &fd_d_model,
      tau, tau_dirs,
      fd_cs, fd_d_cs,
      fd_qdd, fd_qdd_dirs
    );

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);
    // checkConstraintSetsADvsFD(ndirs, ad_cs, ad_d_cs, fd_cs, fd_d_cs);

    CHECK_ARRAY_CLOSE(qdd.data(), ad_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(qdd.data(), fd_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(ad_qdd_dirs.data(), fd_qdd_dirs.data(), nq*nq,
                      array_close_prec);

    q_dirs.setRandom();
    qd_dirs.setRandom();
    tau_dirs.setRandom();
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFForwardDynamicsApplyConstraintForces) {
  srand (421337);
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  ForwardDynamicsApplyConstraintForcesTemplate(*this, 1, 1e-6);
}

// -----------------------------------------------------------------------------

template <typename T>
void ForwardDynamicsContactsKokkevisTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec
) {
  Model   model    = obj.model;
  Model   ad_model = obj.model;
  Model   fd_model = obj.model;

  ADModel ad_d_model = ADModel(ad_model);
  ADModel fd_d_model = ADModel(fd_model);

  ConstraintSet cs = obj.constraint_set;
  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet fd_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = ADConstraintSet(ad_cs, ad_model.dof_count);
  ADConstraintSet fd_d_cs = ADConstraintSet(fd_cs, fd_model.dof_count);

  // set up input quantities
  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;

  VectorNd qdd    = VectorNd::Zero(nq);
  VectorNd ad_qdd = VectorNd::Zero(nq);
  VectorNd fd_qdd = VectorNd::Zero(nq);
  MatrixNd ad_qdd_dirs = MatrixNd::Zero(nq, nq);
  MatrixNd fd_qdd_dirs = MatrixNd::Zero(nq, nq);

  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd qd_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd tau_dirs = MatrixNd::Identity(nq, nq);

  for (unsigned trial = 0; trial < trial_count; trial++) {
    VectorNd q = VectorNd::Random(nq);
    VectorNd qd = VectorNd::Random(nq);
    VectorNd tau = VectorNd::Random(nq);

    ForwardDynamicsContactsKokkevis(model, q, qd, tau, cs, qdd);

    AD::ForwardDynamicsContactsKokkevis(
          ad_model, ad_d_model,
          q, q_dirs,
          qd, qd_dirs,
          tau, tau_dirs,
          ad_cs, ad_d_cs,
          ad_qdd, ad_qdd_dirs);

    FD::ForwardDynamicsContactsKokkevis(
          fd_model, &fd_d_model,
          q, q_dirs,
          qd, qd_dirs,
          tau, tau_dirs,
          fd_cs, fd_d_cs,
          fd_qdd, fd_qdd_dirs);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);
    // checkConstraintSetsADvsFD(ndirs, ad_cs, ad_d_cs, fd_cs, fd_d_cs);

    CHECK_ARRAY_CLOSE(qdd.data(), ad_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(qdd.data(), fd_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(ad_qdd_dirs.data(), fd_qdd_dirs.data(), nq*nq,
                      array_close_prec);

    q_dirs.setRandom();
    qd_dirs.setRandom();
    tau_dirs.setRandom();
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFForwardDynamicsContactsKokkevis) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  ForwardDynamicsContactsKokkevisTemplate(*this, 1, 5e-6);
}


template <typename T>
void ForwardDynamicsContactsKokkevisTemplateFDC(
    T & obj,
    unsigned trial_count,
    double array_close_prec
) {
  Model   model    = obj.model;
  Model   ad_model = obj.model;
  Model   fd_model = obj.model;

  ADModel ad_d_model = ADModel(ad_model);
  ADModel fd_d_model = ADModel(fd_model);

  ConstraintSet cs = obj.constraint_set;
  ConstraintSet ad_cs = obj.constraint_set;
  ConstraintSet fd_cs = obj.constraint_set;
  ADConstraintSet ad_d_cs = ADConstraintSet(ad_cs, ad_model.dof_count);
  ADConstraintSet fd_d_cs = ADConstraintSet(fd_cs, fd_model.dof_count);

  // set up input quantities
  int const nq = ad_model.dof_count;
  unsigned const ndirs = nq;

  VectorNd qdd    = VectorNd::Zero(nq);
  VectorNd ad_qdd = VectorNd::Zero(nq);
  VectorNd fd_qdd = VectorNd::Zero(nq);
  MatrixNd ad_qdd_dirs = MatrixNd::Zero(nq, nq);
  MatrixNd fd_qdd_dirs = MatrixNd::Zero(nq, nq);

  MatrixNd q_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd qd_dirs = MatrixNd::Identity(nq, nq);
  MatrixNd tau_dirs = MatrixNd::Identity(nq, nq);

  for (unsigned trial = 0; trial < trial_count; trial++) {
    VectorNd q = VectorNd::Random(nq);
    VectorNd qd = VectorNd::Random(nq);
    VectorNd tau = VectorNd::Random(nq);

    ForwardDynamicsContactsKokkevis(model, q, qd, tau, cs, qdd);

    AD::ForwardDynamicsContactsKokkevis(
          ad_model, ad_d_model,
          q, q_dirs,
          qd, qd_dirs,
          tau, tau_dirs,
          ad_cs, ad_d_cs,
          ad_qdd, ad_qdd_dirs);

    FDC::ForwardDynamicsContactsKokkevis(
          fd_model, &fd_d_model,
          q, q_dirs,
          qd, qd_dirs,
          tau, tau_dirs,
          fd_cs, fd_d_cs,
          fd_qdd, fd_qdd_dirs);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);
    // checkConstraintSetsADvsFD(ndirs, ad_cs, ad_d_cs, fd_cs, fd_d_cs);

    CHECK_ARRAY_CLOSE(qdd.data(), ad_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(qdd.data(), fd_qdd.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(ad_qdd_dirs.data(), fd_qdd_dirs.data(), nq*nq,
                      array_close_prec);

    q_dirs.setRandom();
    qd_dirs.setRandom();
    tau_dirs.setRandom();
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFForwardDynamicsContactsKokkevis_FDC) {
  // add contacts and bind them to constraint set
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  ForwardDynamicsContactsKokkevisTemplateFDC(*this, 1, 1e-8);
}


// -----------------------------------------------------------------------------

TEST (SolveContactSystemDirectTest) {
  vector<LinearSolver> lss(6);
  lss[0] = LinearSolverPartialPivLU;
  lss[1] = LinearSolverColPivHouseholderQR;
  lss[2] = LinearSolverHouseholderQR;
  lss[3] = LinearSolverPartialPivLU;
  lss[4] = LinearSolverColPivHouseholderQR;
  lss[5] = LinearSolverHouseholderQR;

  for (unsigned i = 0; i < lss.size() * 3; i++) {
    int nc = 3;
    int nm = 11;
    int ndirs = 20;
    LinearSolver ls = lss[i % lss.size()];

    MatrixNd H      = MatrixNd::Identity(nm, nm);
    MatrixNd G      = MatrixNd::Random(nc, nm);
    G.setRandom(nc, nm);
    VectorNd c      = VectorNd::Ones(nm, 1);
    VectorNd gamma  = VectorNd::Ones(nc, 1);
    vector<MatrixNd> H_dirs(ndirs, H);
    vector<MatrixNd> G_dirs(ndirs, G);
    MatrixNd c_dirs = MatrixNd::Random(nm, ndirs);
    MatrixNd gamma_dirs = MatrixNd::Random(nc, ndirs);

    for (int ir = 0; ir < nm; ir++) {
      for (int ic = ir + 1; ic < nm; ic++) {
        H(ir, ic) = MatrixNd::Random(1, 1)(0, 0);
      }
    }

    for (int j = 0; j < nc; j++) {
      H_dirs[j] = MatrixNd::Random(nm, nm);
      G_dirs[j] = MatrixNd::Random(nc, nm);
    }

    MatrixNd A(nm + nc, nm + nc);
    VectorNd b(nm + nc);
    VectorNd x(nm + nc);
    VectorNd qddot(nm);
    VectorNd lambda(nm);
    A.setZero();
    b.setZero();
    SolveConstrainedSystemDirect(H, G, c, gamma, qddot, lambda, A, b, x, ls);

    MatrixNd ad_A(nm + nc, nm + nc);
    VectorNd ad_b(nm + nc);
    VectorNd ad_x(nm + nc);
    ad_A.setZero();
    ad_b.setZero();
    vector<MatrixNd> ad_A_dirs(ndirs, ad_A);
    MatrixNd ad_b_dirs(nm + nc, ndirs);
    MatrixNd ad_x_dirs(nm + nc, ndirs);
    ad_b_dirs.setZero();
    ad_x_dirs.setZero();
    AD::SolveConstrainedSystemDirect(H, H_dirs, G, G_dirs, c, c_dirs, gamma,
                                 gamma_dirs, ad_A, ad_A_dirs, ad_b, ad_b_dirs,
                                 ad_x, ad_x_dirs, ls, ndirs);

    MatrixNd fd_A(nm + nc, nm + nc);
    VectorNd fd_b(nm + nc);
    VectorNd fd_x(nm + nc);
    fd_A.setZero();
    fd_b.setZero();
    vector<MatrixNd> fd_A_dirs(ndirs, ad_A);
    MatrixNd fd_b_dirs(nm + nc, ndirs);
    MatrixNd fd_x_dirs(nm + nc, ndirs);
    fd_b_dirs.setZero();
    fd_x_dirs.setZero();
    FD::SolveConstrainedSystemDirect(H, H_dirs, G, G_dirs, c, c_dirs, gamma,
                                 gamma_dirs, fd_A, fd_A_dirs, fd_b, fd_b_dirs,
                                 fd_x, fd_x_dirs, ls, ndirs);

    CHECK_ARRAY_CLOSE(x.data(), ad_x.data(), nm + nc, 1e-8);
    CHECK_ARRAY_CLOSE(x.data(), fd_x.data(), nm + nc, 1e-8);
    CHECK_ARRAY_CLOSE(b.data(), ad_b.data(), nm + nc, 1e-8);
    CHECK_ARRAY_CLOSE(b.data(), fd_b.data(), nm + nc, 1e-8);
    CHECK_ARRAY_CLOSE(A.data(), fd_A.data(), (nm + nc) * (nm + nc), 1e-8);
    CHECK_ARRAY_CLOSE(A.data(), ad_A.data(), (nm + nc) * (nm + nc), 1e-8);
    for (int idir = 0; idir < ndirs; idir++) {
      CHECK_ARRAY_CLOSE(ad_A_dirs[idir].data(), fd_A_dirs[idir].data(),
                        (nm + nc) * (nm + nc), 1e-4);
    }

    MatrixNd ad_cwise_normalized = ad_x_dirs.cwiseQuotient(ad_x_dirs.cwiseAbs());
    MatrixNd fd_cwise_normalized = fd_x_dirs.cwiseQuotient(fd_x_dirs.cwiseAbs());

    CHECK_ARRAY_CLOSE(ad_cwise_normalized.data(), fd_cwise_normalized.data(),
                      (nm + nc) * ndirs,
                      1e-5);
    CHECK_ARRAY_CLOSE(ad_b_dirs.data(), fd_b_dirs.data(), (nm + nc) * ndirs,
                      1e-4);
  }
}

// -----------------------------------------------------------------------------

template <typename T>
void ComputeContactImpulsesDirectTestTemplate(
    T & obj,
    unsigned trial_count,
    double const array_close_prec) {
  Model   model      = obj.model;
  Model   ad_model   = obj.model;
  Model   fd_model   = obj.model;

  ADModel ad_d_model = ADModel(ad_model);
  ADModel fd_d_model = ADModel(fd_model);

  ConstraintSet   cs      = obj.constraint_set;
  ConstraintSet   ad_cs   = obj.constraint_set;
  ConstraintSet   fd_cs   = obj.constraint_set;
  ADConstraintSet ad_d_cs = ADConstraintSet(ad_cs, ad_model.dof_count);
  ADConstraintSet fd_d_cs = ADConstraintSet(fd_cs, fd_model.dof_count);

  int const nq       = model.dof_count;
  int const ndirs    = 2 * nq;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    VectorNd q         = VectorNd::Random(nq);
    MatrixNd q_dirs    = MatrixNd::Zero(nq, ndirs);
    q_dirs.block(0, 0, nq, nq) = MatrixNd::Identity(nq, nq);

    VectorNd qd        = VectorNd::Random(nq);
    MatrixNd qd_dirs   = MatrixNd::Zero(nq, ndirs);
    qd_dirs.block(0, nq, nq, nq) = MatrixNd::Identity(nq, nq);

    VectorNd qd_plus(nq);
    VectorNd ad_qd_plus(nq);
    MatrixNd ad_qd_plus_dirs(nq, ndirs);
    VectorNd fd_qd_plus(nq);
    MatrixNd fd_qd_plus_dirs(nq, ndirs);

    ComputeConstraintImpulsesDirect (model, q, qd, cs, qd_plus);
    AD::ComputeConstraintImpulsesDirect(ad_model, ad_d_model,
                                     q, q_dirs, qd, qd_dirs,
                                     ad_cs, ad_d_cs,
                                     ad_qd_plus, ad_qd_plus_dirs);
    FD::ComputeConstraintImpulsesDirect (fd_model, &fd_d_model,
                                         q, q_dirs, qd, qd_dirs,
                                         fd_cs, &fd_d_cs,
                                         fd_qd_plus, fd_qd_plus_dirs);

    checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);
    // checkConstraintSetsADvsFD(ndirs, ad_cs, ad_d_cs, fd_cs, fd_d_cs);

    CHECK_ARRAY_CLOSE(qd_plus.data(), ad_qd_plus.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(qd_plus.data(), fd_qd_plus.data(), nq, array_close_prec);
    CHECK_ARRAY_CLOSE(ad_qd_plus_dirs.data(), fd_qd_plus_dirs.data(),
                      nq * ndirs, array_close_prec);
  }
}

TEST_FIXTURE (FixedBase6DoF, FixedBase6DoFComputeContactImpulsesDirectTest) {
  constraint_set.AddContactConstraint (
        contact_body_id, Vector3d (1., 0., 0.), contact_normal);
  constraint_set.AddContactConstraint (
        contact_body_id, Vector3d (0., 1., 0.), contact_normal);
  constraint_set.Bind (model);
  ComputeContactImpulsesDirectTestTemplate(*this, 5, 1e-5);
}

// -----------------------------------------------------------------------------

