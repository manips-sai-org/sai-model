#include <UnitTest++.h>

#include <iomanip>
#include <iostream>
#include <random>

#include "Fixtures.h"
#include "rbdl_utilsFD.h"
#include "rbdl_utilsAD.h"

#include "SpatialAlgebraOperatorsAD.h"
#include "ModelCheckADvsFD.h"
#include "Human36Fixture.h"

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;
using namespace RigidBodyDynamics::Utils;

const double TEST_PREC = 1.0e-12;

using namespace std;

// -----------------------------------------------------------------------------

SpatialRigidBodyInertia operator* (
    double s,
    SpatialRigidBodyInertia const & rhs) {
  SpatialRigidBodyInertia res = rhs;

  res.m *= s;
  res.h *= s;
  res.Ixx *= s;
  res.Iyx *= s;
  res.Iyy *= s;
  res.Izx *= s;
  res.Izy *= s;
  res.Izz *= s;

  return res;
}

SpatialRigidBodyInertia operator* (
    SpatialRigidBodyInertia const & lhs,
    double s) {
  return s * lhs;
}

SpatialRigidBodyInertia operator/ (
    SpatialRigidBodyInertia const & lhs,
    double s) {
  return (1. / s) * lhs;
}

SpatialRigidBodyInertia operator+ (
    SpatialRigidBodyInertia const & lhs,
    SpatialRigidBodyInertia const & rhs) {
  SpatialRigidBodyInertia res;

  res.m   = lhs.m   + rhs.m;
  res.h   = lhs.h   + rhs.h;
  res.Ixx = lhs.Ixx + rhs.Ixx;
  res.Iyx = lhs.Iyx + rhs.Iyx;
  res.Iyy = lhs.Iyy + rhs.Iyy;
  res.Izx = lhs.Izx + rhs.Izx;
  res.Izy = lhs.Izy + rhs.Izy;
  res.Izz = lhs.Izz + rhs.Izz;

  return res;
}

SpatialRigidBodyInertia operator -(
    SpatialRigidBodyInertia const & lhs,
    SpatialRigidBodyInertia const & rhs) {
  SpatialRigidBodyInertia res = lhs;

  res.m   -= rhs.m;
  res.h   -= rhs.h;
  res.Ixx -= rhs.Ixx;
  res.Iyx -= rhs.Iyx;
  res.Iyy -= rhs.Iyy;
  res.Izx -= rhs.Izx;
  res.Izy -= rhs.Izy;
  res.Izz -= rhs.Izz;

  return res;
}


void FDapplyTransposeSTSI (
    SpatialTransform const & st,
    vector<SpatialTransform> const & st_dirs,
    SpatialRigidBodyInertia const & srbi,
    vector<SpatialRigidBodyInertia> const & srbi_dirs,
    SpatialRigidBodyInertia & res,
    vector<SpatialRigidBodyInertia> & res_dirs) {

  unsigned const ndirs = st_dirs.size();

  assert(ndirs == srbi_dirs.size());
  assert(ndirs == res_dirs.size());

  // nominal code
  res = st.applyTranspose(srbi);

  double h = 1e-8;

  // derivative code
  for (unsigned i = 0; i < ndirs; i++) {

    SpatialTransform st_h;
    st_h.r    = st.r + h * st_dirs[i].r;
    st_h.E    = st.E + h * st_dirs[i].E;

    SpatialRigidBodyInertia srbi_h = srbi + h * srbi_dirs[i];
    res_dirs[i]  = st_h.applyTranspose(srbi_h);
    res_dirs[i]  = res_dirs[i] - res;
    res_dirs[i]  = res_dirs[i] / h;
  }
}

TEST(CheckApplyTranspose) {
  std::uniform_real_distribution<> dist(0.0, 1.0);
  std::default_random_engine gen;

  for (int trial = 0; trial < 100; trial++){
    Matrix3d E = Matrix3d::Random();
    Vector3d r = Vector3d::Random();

    double m   = dist(gen);
    Vector3d h = m * Vector3d::Random();
    double Ixx = dist(gen);
    double Ixy = dist(gen);
    double Ixz = dist(gen);
    double Iyy = dist(gen);
    double Iyz = dist(gen);
    double Izz = dist(gen);

    Matrix3d I (Ixx, Ixy, Ixz,
                Ixy, Iyy, Iyz,
                Ixz, Iyz, Izz);
    SpatialTransform        st(E, r);
    SpatialRigidBodyInertia srbi(m, h, I);

    unsigned const ndirs = 22;

    vector<SpatialTransform> st_dirs(ndirs);
    vector<SpatialRigidBodyInertia> srbi_dirs(ndirs);
    for (unsigned i = 0; i < 9; i++) {
      Matrix3d E = Matrix3d::Zero();
      E(i / 3, i % 3) = 1.0;
      st_dirs[i] = SpatialTransform(E, Vector3d::Zero());
    }
    for (unsigned i = 0; i < 3; i++) {
      Vector3d r = Vector3d::Zero();
      r(i) = 1.0;
      st_dirs[9 + i] = SpatialTransform(Matrix3d::Zero(), r);
    }
    srbi_dirs[12] = SpatialRigidBodyInertia(1.0, Vector3d::Zero(), Matrix3d::Zero());
    for (unsigned i = 0; i < 3; i++) {
      Vector3d h = Vector3d::Zero();
      h(i) = 1.0;
      srbi_dirs[13 + i] = SpatialRigidBodyInertia(0.0, h, Matrix3d::Zero());
    }
    srbi_dirs[16] = SpatialRigidBodyInertia(0.0, Vector3d::Zero(), Matrix3d(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    srbi_dirs[17] = SpatialRigidBodyInertia(0.0, Vector3d::Zero(), Matrix3d(0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    srbi_dirs[18] = SpatialRigidBodyInertia(0.0, Vector3d::Zero(), Matrix3d(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0));
    srbi_dirs[19] = SpatialRigidBodyInertia(0.0, Vector3d::Zero(), Matrix3d(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0));
    srbi_dirs[20] = SpatialRigidBodyInertia(0.0, Vector3d::Zero(), Matrix3d(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0));
    srbi_dirs[21] = SpatialRigidBodyInertia(0.0, Vector3d::Zero(), Matrix3d(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0));

    SpatialRigidBodyInertia ad_res;
    vector<SpatialRigidBodyInertia> ad_res_dirs(ndirs);
    AD::applyTransposeSTSI(ndirs, st, st_dirs, srbi, srbi_dirs, ad_res, ad_res_dirs);

    SpatialRigidBodyInertia fd_res;
    vector<SpatialRigidBodyInertia> fd_res_dirs(ndirs);
    FDapplyTransposeSTSI(st, st_dirs, srbi, srbi_dirs, fd_res, fd_res_dirs);

    CHECK_EQUAL(ad_res.m, fd_res.m);
    CHECK_EQUAL(ad_res.h, fd_res.h);
    CHECK_EQUAL(ad_res.Ixx, fd_res.Ixx);
    CHECK_EQUAL(ad_res.Iyx, fd_res.Iyx);
    CHECK_EQUAL(ad_res.Iyy, fd_res.Iyy);
    CHECK_EQUAL(ad_res.Izx, fd_res.Izx);
    CHECK_EQUAL(ad_res.Izy, fd_res.Izy);
    CHECK_EQUAL(ad_res.Izz, fd_res.Izz);

    for (unsigned j = 0; j < ndirs; j++) {
      CHECK_CLOSE(ad_res_dirs[j].m, fd_res_dirs[j].m, 1e-6);
      CHECK_ARRAY_CLOSE(ad_res_dirs[j].h.data(), fd_res_dirs[j].h.data(), 3, 1e-6);
      CHECK_CLOSE(ad_res_dirs[j].Ixx, fd_res_dirs[j].Ixx, 1e-6);
      CHECK_CLOSE(ad_res_dirs[j].Iyx, fd_res_dirs[j].Iyx, 1e-6);
      CHECK_CLOSE(ad_res_dirs[j].Iyy, fd_res_dirs[j].Iyy, 1e-6);
      CHECK_CLOSE(ad_res_dirs[j].Izx, fd_res_dirs[j].Izx, 1e-6);
      CHECK_CLOSE(ad_res_dirs[j].Izy, fd_res_dirs[j].Izy, 1e-6);
      CHECK_CLOSE(ad_res_dirs[j].Izz, fd_res_dirs[j].Izz, 1e-6);
    }
  }
}

// -----------------------------------------------------------------------------

inline void FDapplySTSV(
    unsigned ndirs,
    SpatialTransform const & st,
    vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    vector<SpatialVector> const & sv_dirs,
    SpatialVector & res,
    vector<SpatialVector> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.size());
  assert(ndirs <= res_dirs.size());

  res = st.apply(sv);
  double h = 1e-8;

  for (unsigned idir = 0; idir < ndirs; idir++) {
    SpatialTransform sth;
    sth.E = st.E + h * st_dirs[idir].E;
    sth.r = st.r + h * st_dirs[idir].r;
    SpatialVector svh;
    svh = sv + h * sv_dirs[idir];
    res_dirs[idir] = (sth.apply(svh) - res) / h;
  }
}

TEST(CheckApplySTSV) {
  uniform_real_distribution<> dist(0.0, 1.0);
  default_random_engine gen;

  for (int i = 0; i < 100; i++){
    unsigned ndirs = 18;
    Matrix3d E = Matrix3d::Random();
    Vector3d r = Vector3d::Random();
    SpatialTransform st(E, r);
    SpatialVector sv = SpatialVector::Random();
    vector<SpatialTransform> st_dirs(ndirs, SpatialTransform::Zero());
    vector<SpatialVector> sv_dirs(ndirs, SpatialVector::Zero());

    SpatialVector ad_res = SpatialVector::Random();
    vector<SpatialVector> ad_res_dirs(ndirs, SpatialVector::Random());

    SpatialVector fd_res = ad_res;
    vector<SpatialVector> fd_res_dirs = ad_res_dirs;

    unsigned idir = 0;
    for (; idir < 9u; idir++) {
      st_dirs[idir].E(idir % 3, idir / 3) = 1.0;
    }
    for (; idir < 12u; idir++) {
      st_dirs[idir].r(idir - 9u) = 1.0;
    }
    for (; idir < ndirs; idir++) {
      sv_dirs[idir](idir - 12u) = 1.0;
    }

    AD::applySTSV(
          ndirs,
          st, st_dirs,
          sv, sv_dirs,
          ad_res, ad_res_dirs);

    FDapplySTSV(
          ndirs,
          st, st_dirs,
          sv, sv_dirs,
          fd_res, fd_res_dirs);

    CHECK_ARRAY_EQUAL(ad_res.data(), fd_res.data(), 6);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      CHECK_ARRAY_CLOSE(ad_res_dirs[idir].data(), fd_res_dirs[idir].data(), 6, 1e-7);
    }
  }
}

// -----------------------------------------------------------------------------

inline void FDapplyTransposeSTSV(
    unsigned ndirs,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    std::vector<SpatialVector> const & sv_dirs,
    SpatialVector & res,
    std::vector<SpatialVector> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.size());
  assert(ndirs <= res_dirs.size());

  res = st.applyTranspose(sv);
  double h = 1e-8;

  for (unsigned idir = 0; idir < ndirs; idir++) {
    SpatialTransform sth;
    sth.E = st.E + h * st_dirs[idir].E;
    sth.r = st.r + h * st_dirs[idir].r;
    SpatialVector svh;
    svh = sv + h * sv_dirs[idir];
    res_dirs[idir] = (sth.applyTranspose(svh) - res) / h;
  }
}

TEST(CheckApplyTransposeSTSV) {
  for (int i = 0; i < 100; i++){
    unsigned ndirs = 18;
    Matrix3d E = Matrix3d::Random();
    Vector3d r = Vector3d::Random();
    SpatialTransform st(E, r);
    SpatialVector sv = SpatialVector::Random();
    vector<SpatialTransform> st_dirs(ndirs, SpatialTransform::Zero());
    vector<SpatialVector> sm_dirs(ndirs, SpatialVector::Zero());

    SpatialVector ad_res = SpatialVector::Zero();
    vector<SpatialVector> ad_res_dirs(ndirs, SpatialVector::Zero());

    SpatialVector fd_res = SpatialVector::Zero();
    vector<SpatialVector> fd_res_dirs(ndirs, SpatialVector::Zero());

    unsigned idir = 0;
    for (; idir < 9u; idir++) {
      st_dirs[idir].E(idir % 3, idir / 3) = 1.0;
    }
    for (; idir < 12u; idir++) {
      st_dirs[idir].r(idir - 9u) = 1.0;
    }
    for (; idir < ndirs; idir++) {
      sm_dirs[idir](idir - 12u) = 1.0;
    }
    AD::applyTransposeSTSV(
          ndirs,
          st, st_dirs,
          sv, sm_dirs,
          ad_res, ad_res_dirs);

    FDapplyTransposeSTSV(
          ndirs,
          st, st_dirs,
          sv, sm_dirs,
          fd_res, fd_res_dirs);

    CHECK_ARRAY_CLOSE(ad_res.data(), fd_res.data(), 6, 1e-7);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      CHECK_ARRAY_CLOSE(ad_res_dirs[idir].data(), fd_res_dirs[idir].data(), 6, 1e-7);
    }
  }
}

// -----------------------------------------------------------------------------

inline void FDaddApplyAdjointSTSV (
    unsigned ndirs,
    double dscal,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    std::vector<SpatialVector> const & sv_dirs,
    SpatialVector & res,
    std::vector<SpatialVector> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.size());
  assert(ndirs <= res_dirs.size());

  SpatialVector offset = res;
  res += dscal * st.toMatrixAdjoint() * sv;
  double h = 1e-8;
  for (unsigned idir = 0; idir < ndirs; idir++) {
    SpatialTransform sth = st;
    sth.E = st.E + h * st_dirs[idir].E;
    sth.r = st.r + h * st_dirs[idir].r;
    SpatialVector svh;
    svh = sv + h * sv_dirs[idir];
    res_dirs[idir] += (dscal * (sth.toMatrixAdjoint() * svh) - (res - offset)) / h;
  }
}

TEST(CheckAddApplyAdjointSTSV) {
  std::uniform_real_distribution<> dist(-2.0, 2.0);
  std::default_random_engine gen;
  for (int i = 0; i < 100; i++){
    unsigned ndirs = 18;
    Matrix3d E = Matrix3d::Random();
    Vector3d r = Vector3d::Random();
    SpatialTransform st(E, r);
    SpatialVector sv = SpatialVector::Random();
    vector<SpatialTransform> st_dirs(ndirs, SpatialTransform::Zero());
    vector<SpatialVector> sv_dirs(ndirs, SpatialVector::Zero());

    SpatialVector ad_res = SpatialVector::Random();
    vector<SpatialVector> ad_res_dirs(ndirs);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_res_dirs[idir].setRandom();
    }

    SpatialVector fd_res = ad_res;
    vector<SpatialVector> fd_res_dirs = ad_res_dirs;

    unsigned idir = 0;
    for (; idir < 9u; idir++) {
      st_dirs[idir].E(idir % 3, idir / 3) = 1.0;
    }
    for (; idir < 12u; idir++) {
      st_dirs[idir].r(idir - 9u) = 1.0;
    }
    for (; idir < ndirs; idir++) {
      sv_dirs[idir](idir - 12u) = 1.0;
    }

    double dscal = dist(gen);

    AD::addApplyAdjointSTSV(
          ndirs,
          dscal,
          st, st_dirs,
          sv, sv_dirs,
          ad_res, ad_res_dirs);

    FDaddApplyAdjointSTSV(
          ndirs,
          dscal,
          st, st_dirs,
          sv, sv_dirs,
          fd_res, fd_res_dirs);

    CHECK_ARRAY_CLOSE(ad_res.data(), fd_res.data(), 6, 1e-7);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      CHECK_ARRAY_CLOSE(ad_res_dirs[idir].data(), fd_res_dirs[idir].data(), 6, 1e-6);
    }
  }
}

// -----------------------------------------------------------------------------

inline void FDaddSqrFormSTSM_noalias(
    unsigned ndirs,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialMatrix const & sm,
    std::vector<SpatialMatrix> const & sm_dirs,
    SpatialMatrix & res,
    std::vector<SpatialMatrix> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sm_dirs.size());
  assert(ndirs <= res_dirs.size());

  SpatialMatrix offset = res;
  res.noalias() += st.toMatrixTranspose() * sm * st.toMatrix() ;
  double h = 1e-8;

  for (unsigned idir = 0; idir < ndirs; idir++) {
    SpatialTransform sth = st;
    sth.E = st.E + h * st_dirs[idir].E;
    sth.r = st.r + h * st_dirs[idir].r;
    SpatialMatrix smh;
    smh = sm + h * sm_dirs[idir];

    res_dirs[idir].noalias() += ((sth.toMatrixTranspose() * smh * sth.toMatrix()) - (res - offset)) / h;
  }
}

TEST(CheckAddSqrFormSTSM_noalias) {
  for (unsigned trial = 0; trial < 100; trial++){
    unsigned ndirs = 48;
    Matrix3d E = Matrix3d::Random();
    Vector3d r = Vector3d::Random();
    SpatialTransform st(E, r);
    SpatialMatrix sm = SpatialMatrix::Random();
    vector<SpatialTransform> st_dirs(ndirs, SpatialTransform::Zero());
    vector<SpatialMatrix> sm_dirs(ndirs, SpatialMatrix::Zero());

    SpatialMatrix ad_res = SpatialMatrix::Random();
    vector<SpatialMatrix> ad_res_dirs(ndirs, SpatialMatrix::Zero());

    SpatialMatrix fd_res = ad_res;
    vector<SpatialMatrix> fd_res_dirs(ndirs, SpatialMatrix::Zero());

    unsigned idir = 0;
    for (; idir < 9u; idir++) {
      st_dirs[idir].E(idir % 3, idir / 3) = 1.0;
    }
    for (; idir < 12u; idir++) {
      st_dirs[idir].r(idir - 9u) = 1.0;
    }
    for (; idir < ndirs; idir++) {
      sm_dirs[idir]((idir - 12u) / 6, (idir - 12u) % 6) = 1.0;
    }

    AD::addSqrFormSTSM_noalias(
          ndirs,
          st, st_dirs,
          sm, sm_dirs,
          ad_res, ad_res_dirs);

    FDaddSqrFormSTSM_noalias(
          ndirs,
          st, st_dirs,
          sm, sm_dirs,
          fd_res, fd_res_dirs);

    CHECK_ARRAY_CLOSE(ad_res.data(), fd_res.data(), 6, 1e-7);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      CHECK_ARRAY_CLOSE(ad_res_dirs[idir].data(), fd_res_dirs[idir].data(), 36, 1e-6);
    }
  }
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

template <typename T>
void CalcPotentialEnergyADTestTemplate(T & obj) {
  Model & model = obj.model;
  ADModel ad_model = ADModel(model);

  VectorNd & q = obj.q;
  VectorNd & qdot = obj.qdot;

  q[0] = 0.3;
  q[1] = -0.2;
  qdot[0] = .1;
  qdot[1] = -.15;

  int nrows = model.dof_count;
  int ndirs = 2 * model.dof_count;

  MatrixNd q_dirs = MatrixNd::Zero(nrows, ndirs);
  MatrixNd qdot_dirs = MatrixNd::Zero(nrows, ndirs);

  q_dirs.block(0, 0, nrows, model.dof_count)
      = MatrixNd::Identity(nrows, model.dof_count);
  q_dirs.block(0, model.dof_count, nrows, model.dof_count)
      = MatrixNd::Zero(nrows, model.dof_count);
  qdot_dirs.block(0, 0, nrows, model.dof_count)
      = MatrixNd::Zero(nrows, model.dof_count);
  qdot_dirs.block(0, model.dof_count, nrows, model.dof_count)
      = MatrixNd::Identity(nrows, model.dof_count);

  MatrixNd fd_pote = MatrixNd::Zero(1, ndirs);
  FD::CalcPotentialEnergy(model, q, q_dirs, fd_pote);

  MatrixNd ad_pote = MatrixNd::Zero(1, ndirs);
  AD::CalcPotentialEnergy(model, ad_model, q, q_dirs, ad_pote, true);

  CHECK_ARRAY_CLOSE(fd_pote.data(), ad_pote.data(), 1 * ndirs, 1e-6);
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcPotentialEnergyADTest) {
  CalcPotentialEnergyADTestTemplate(*this);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcPotentialEnergyADTest) {
  CalcPotentialEnergyADTestTemplate(*this);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcPotentialEnergyADTest) {
  CalcPotentialEnergyADTestTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcPotentialEnergyADTest) {
  CalcPotentialEnergyADTestTemplate(*this);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcPotentialEnergyADTest) {
  CalcPotentialEnergyADTestTemplate(*this);
}

// -----------------------------------------------------------------------------

template <typename T>
void CalcKineticEnergyADTestTemplate(
    T & obj,
    unsigned trial_count,
    double array_close_prec) {
  Model & model = obj.model;
  ADModel ad_model = ADModel(model);

  VectorNd & q = obj.q;
  VectorNd & qdot = obj.qdot;

  int nrows = model.dof_count;
  int ndirs = 2 * model.dof_count;

  for (unsigned trial = 0; trial < trial_count; trial++) {
    q.setRandom();
    qdot.setRandom();
    MatrixNd q_dirs = MatrixNd::Random(nrows, ndirs);
    MatrixNd qdot_dirs = MatrixNd::Random(nrows, ndirs);

    MatrixNd fd_kine = MatrixNd::Zero(1, ndirs);
    Utils::FD::CalcKineticEnergy(model, q, q_dirs, qdot, qdot_dirs, fd_kine);

    MatrixNd ad_kine = MatrixNd::Zero(1, ndirs);
    AD::CalcKineticEnergy(model, ad_model, q, q_dirs, qdot, qdot_dirs, ad_kine, true);

    CHECK_ARRAY_CLOSE(fd_kine.data(), ad_kine.data(), 1 * ndirs, array_close_prec);
  }
}

TEST_FIXTURE ( CartPendulum, CartPendulumCalcKineticEnergyADTest) {
  CalcKineticEnergyADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm2DofX, Arm2DofXCalcKineticEnergyADTest) {
  CalcKineticEnergyADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm2DofZ, Arm2DofZCalcKineticEnergyADTest) {
  CalcKineticEnergyADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm3DofXZYp, Arm3DofXZYpCalcKineticEnergyADTest) {
  CalcKineticEnergyADTestTemplate(*this, 10, 1e-5);
}

TEST_FIXTURE ( Arm3DofXZZp, Arm3DofXZZpCalcKineticEnergyADTest) {
  CalcKineticEnergyADTestTemplate(*this, 10, 1e-5);
}

// -----------------------------------------------------------------------------

template<typename T>
void CalcCenterOfMassADTestTemplate(T & obj) {
  Model ad_model = obj.model;
  Model fd_model = obj.model;

  ADModel ad_d_model = ADModel(ad_model);
  ADModel fd_d_model = ADModel(fd_model);

  VectorNd & q = obj.q;
  VectorNd & qdot = obj.qdot;

  q[0]    =  M_PI / 4.;
  q[1]    = -M_PI / 3.;
  qdot[0] = -.2;
  qdot[1] =  .11;

  int nrows = obj.model.dof_count;
  int ndirs = 2 * obj.model.dof_count;

  MatrixNd q_dirs = MatrixNd::Zero(nrows, ndirs);
  MatrixNd qdot_dirs = MatrixNd::Zero(nrows, ndirs);

  q_dirs.block(0, 0, nrows, obj.model.dof_count)
      = MatrixNd::Identity(nrows, obj.model.dof_count);
  q_dirs.block(0, obj.model.dof_count, nrows, obj.model.dof_count)
      = MatrixNd::Zero(nrows, obj.model.dof_count);
  qdot_dirs.block(0, 0, nrows, obj.model.dof_count)
      = MatrixNd::Zero(nrows, obj.model.dof_count);
  qdot_dirs.block(0, obj.model.dof_count, nrows, obj.model.dof_count)
      = MatrixNd::Identity(nrows, obj.model.dof_count);

  double   ad_mass;
  Vector3d ad_com;
  MatrixNd ad_d_com(3, q_dirs.cols());
  Vector3d ad_comVelocity(0., 0., 0.);
  MatrixNd ad_d_comVelocity(3, ndirs);
  Vector3d ad_angMomentum(0., 0., 0.);
  MatrixNd ad_d_angMomentum(3, ndirs);
  AD::CalcCenterOfMass(ad_model, ad_d_model, q, q_dirs, qdot, qdot_dirs, ad_mass,
                       ad_com, ad_d_com,
                       &ad_comVelocity, &ad_d_comVelocity,
                       &ad_angMomentum, &ad_d_angMomentum, true);

  double   fd_mass;
  Vector3d fd_com;
  MatrixNd fd_d_com(3, ndirs);
  Vector3d fd_comVelocity(0., 0., 0.);
  MatrixNd fd_d_comVelocity(3, ndirs);
  Vector3d fd_angMomentum(0., 0., 0.);
  MatrixNd fd_d_angMomentum(3, ndirs);
  FD::CalcCenterOfMass(fd_model, &fd_d_model, q, q_dirs, qdot, qdot_dirs, fd_mass, fd_com,
                       fd_d_com, &fd_comVelocity, &fd_d_comVelocity,
                       &fd_angMomentum, &fd_d_angMomentum);

  checkModelsADvsFD(ndirs, ad_model, ad_d_model, fd_model, fd_d_model);

  CHECK_EQUAL(ad_mass, fd_mass);
  CHECK_ARRAY_CLOSE(ad_com.data(), fd_com.data(), 3, TEST_PREC);
  CHECK_ARRAY_CLOSE(ad_d_com.data(), fd_d_com.data(), 3 * ndirs, 1e-8);
  CHECK_ARRAY_CLOSE(ad_comVelocity.data(), fd_comVelocity.data(), 3, 1e-8);
  CHECK_ARRAY_CLOSE(ad_d_comVelocity.data(), fd_d_comVelocity.data(), 3 * ndirs, 1e-8);
  CHECK_ARRAY_CLOSE(ad_angMomentum.data(), fd_angMomentum.data(), 3, 1e-8);
  CHECK_ARRAY_CLOSE(ad_d_angMomentum.data(), fd_d_angMomentum.data(), 3 * ndirs, 1e-8);
}

TEST_FIXTURE( CartPendulum, CartPendulumCalcCenterOfMass) {
  CalcCenterOfMassADTestTemplate(*this);
}

TEST_FIXTURE( Arm2DofX, Arm2DofXCalcCenterOfMass) {
  CalcCenterOfMassADTestTemplate(*this);
}

TEST_FIXTURE( Arm2DofZ, Arm2DofZCalcCenterOfMass) {
  CalcCenterOfMassADTestTemplate(*this);
}

TEST_FIXTURE( Arm3DofXZYp, Arm3DofXZYpCalcCenterOfMass) {
  CalcCenterOfMassADTestTemplate(*this);
}

TEST_FIXTURE( Arm3DofXZZp, Arm3DofXZZpCalcCenterOfMass) {
  CalcCenterOfMassADTestTemplate(*this);
}

// -----------------------------------------------------------------------------
