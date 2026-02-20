#include <UnitTest++.h>

#include <iostream>

#include "rbdl/rbdl_mathutils.h"

#include "rbdl_mathutilsAD.h"
#include "rbdl_mathutilsFD.h"

#include "SpatialAlgebraOperatorsAD.h"

#include "rbdl/Logging.h"

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const double TEST_PREC = 1.0e-12;

// -----------------------------------------------------------------------------

TEST (crossm_v1v2_ADvsAnalyticTest) {
  unsigned ndirs = 20;

  SpatialVector v1 = SpatialVector::Random();
  SpatialVector v2 = SpatialVector::Random();
  MatrixNd v1_dirs = MatrixNd::Random(6, ndirs);
  MatrixNd v2_dirs = MatrixNd::Random(6, ndirs);

  vector<SpatialVector> ad_res (ndirs, SpatialVector::Zero());
  vector<SpatialVector> an_res (ndirs, SpatialVector::Zero());

  for (unsigned idir = 0; idir < ndirs; idir++) {
    SpatialVector v1_dir = v1_dirs.block(0, idir, 6, 1);
    SpatialVector v2_dir = v2_dirs.block(0, idir, 6, 1);

    ad_res[idir] = AD::crossm (v1, v1_dir, v2, v2_dir);
    an_res[idir] = crossm (v1_dir, v2) + crossm (v1, v2_dir);
    CHECK_ARRAY_CLOSE (
          an_res[idir].data(), ad_res[idir].data(), 6, TEST_PREC);
  }
}

TEST (crossm_v1v2_ADvsFDTest) {
  double TEST_PREC = 1e-07;
  unsigned int ndirs = 20;

  SpatialVector v1 = SpatialVector::Random();
  SpatialVector v2 = SpatialVector::Random();
  MatrixNd v1_dirs = MatrixNd::Random(6, ndirs);
  MatrixNd v2_dirs = MatrixNd::Random(6, ndirs);

  vector<SpatialVector> ad_res (ndirs, SpatialVector::Zero());
  vector<SpatialVector> fd_res (ndirs, SpatialVector::Zero());

  for (unsigned int idir = 0; idir < ndirs; idir++) {
    SpatialVector v1_dir = v1_dirs.block(0, idir, 6, 1);
    SpatialVector v2_dir = v2_dirs.block(0, idir, 6, 1);

    ad_res[idir] = AD::crossm (v1, v1_dir, v2, v2_dir);
    fd_res[idir] = FD::crossm (v1, v1_dir, v2, v2_dir);
    CHECK_ARRAY_CLOSE (
          fd_res[idir].data(), ad_res[idir].data(), 6, TEST_PREC);
  }
}

TEST (crossm_v_ADvsAnalyticTest) {
  double TEST_PREC = 1e-08;
  unsigned ndirs = 20;

  MatrixNd v_dirs = MatrixNd::Random(6, ndirs);
  vector<SpatialMatrix> ad_res (ndirs, SpatialMatrix::Zero());
  vector<SpatialMatrix> an_res (ndirs, SpatialMatrix::Zero());
  for (unsigned idir = 0; idir < ndirs; idir++) {
    SpatialVector v_dir = v_dirs.block(0, idir, 6, 1);
    ad_res[idir] = AD::crossm (v_dir);
    an_res[idir] = crossm (v_dir);
    CHECK_ARRAY_CLOSE (an_res[idir].data(), ad_res[idir].data(), 6, TEST_PREC);
  }
}

TEST (crossm_v_ADvsFDTest) {
  double TEST_PREC = 1e-08;
  unsigned int ndirs = 20;

  SpatialVector v = SpatialVector::Random();
  MatrixNd v_dirs = MatrixNd::Random(6, ndirs);
  vector<SpatialMatrix> ad_res (ndirs, SpatialMatrix::Zero());
  vector<SpatialMatrix> fd_res (ndirs, SpatialMatrix::Zero());

  for (unsigned idir = 0; idir < ndirs; idir++) {
    SpatialVector v_dir = v_dirs.block(0, idir, 6, 1);
    ad_res[idir] = AD::crossm (v_dir);
    fd_res[idir] = FD::crossm (v, v_dir);
    CHECK_ARRAY_CLOSE (fd_res[idir].data(), ad_res[idir].data(), 6, TEST_PREC);
  }
}


TEST (crossf_v1v2_ADvsFDTest) {
  double TEST_PREC = 1e-07;
  unsigned int ndirs = 20;

  SpatialVector v1 = SpatialVector::Random();
  SpatialVector v2 = SpatialVector::Random();
  MatrixNd v1_dirs = MatrixNd::Random(6, ndirs);
  MatrixNd v2_dirs = MatrixNd::Random(6, ndirs);

  vector<SpatialVector> ad_res (ndirs, SpatialVector::Zero());
  vector<SpatialVector> fd_res (ndirs, SpatialVector::Zero());

  for (unsigned int idir = 0; idir < ndirs; idir++) {
    SpatialVector v1_dir = v1_dirs.block(0, idir, 6, 1);
    SpatialVector v2_dir = v2_dirs.block(0, idir, 6, 1);

    ad_res[idir] = AD::crossf (v1, v1_dir, v2, v2_dir);
    fd_res[idir] = FD::crossf (v1, v1_dir, v2, v2_dir);
    CHECK_ARRAY_CLOSE (
          fd_res[idir].data(), ad_res[idir].data(), 6, TEST_PREC);
  }
}


TEST (addApplyAdjointSTSV_vecvec_ADvsFDTest) {
  double TEST_PREC = 1e-6;
  unsigned int ndirs = 20;

  SpatialTransform st(Matrix3d::Random(), Vector3d::Random());
  vector<SpatialTransform> st_dirs(ndirs);
  SpatialVector sv = SpatialVector::Random();
  vector<SpatialVector> sv_dirs(ndirs);
  SpatialVector res_ad = SpatialVector::Random();
  vector<SpatialVector> res_ad_dirs(ndirs);
  SpatialVector res_fd = res_ad;
  vector<SpatialVector> res_fd_dirs(ndirs);
  for (unsigned idir = 0; idir < ndirs; idir++) {
    st_dirs[idir] = SpatialTransform(Matrix3d::Random(), Vector3d::Random());
    sv_dirs[idir] = SpatialVector::Random();
    res_ad_dirs[idir] = SpatialVector::Random();
    res_fd_dirs[idir] = res_ad_dirs[idir];
  }

  double dscal = -3.1415;
  AD::addApplyAdjointSTSV(ndirs, dscal, st, st_dirs, sv, sv_dirs, res_ad, res_ad_dirs);
  FD::addApplyAdjointSTSV(ndirs, dscal, st, st_dirs, sv, sv_dirs, res_fd, res_fd_dirs);

  CHECK_ARRAY_CLOSE (res_ad, res_fd, 6, TEST_PREC);
  for (unsigned int idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE (
          res_ad_dirs[idir].data(), res_fd_dirs[idir].data(), 6, TEST_PREC);
  }
}


TEST (addApplyAdjointSTSV_mat_ADvsFDTest) {
  double TEST_PREC = 1e-6;
  unsigned int ndirs = 20;

  SpatialTransform st(Matrix3d::Random(), Vector3d::Random());
  vector<SpatialTransform> st_dirs(ndirs);
  SpatialVector sv = SpatialVector::Random();
  MatrixNd sv_dirs = MatrixNd::Random(6, ndirs);
  SpatialVector res_ad = SpatialVector::Random();
  vector<SpatialVector> res_ad_dirs(ndirs);
  SpatialVector res_fd = res_ad;
  vector<SpatialVector> res_fd_dirs(ndirs);
  for (unsigned idir = 0; idir < ndirs; idir++) {
    st_dirs[idir] = SpatialTransform(Matrix3d::Random(), Vector3d::Random());
    res_ad_dirs[idir] = SpatialVector::Random();
    res_fd_dirs[idir] = res_ad_dirs[idir];
  }

  double dscal = -3.1415;
  AD::addApplyAdjointSTSV(ndirs, dscal, st, st_dirs, sv, sv_dirs, res_ad, res_ad_dirs);
  FD::addApplyAdjointSTSV(ndirs, dscal, st, st_dirs, sv, sv_dirs, res_fd, res_fd_dirs);

  CHECK_ARRAY_CLOSE (res_ad, res_fd, 6, TEST_PREC);
  for (unsigned int idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE (
          res_ad_dirs[idir].data(), res_fd_dirs[idir].data(), 6, TEST_PREC);
  }
}


// -----------------------------------------------------------------------------

void FDinverse(
    unsigned ndirs,
    SpatialTransform const & st,
    vector<SpatialTransform> const & st_dirs,
    SpatialTransform & res,
    vector<SpatialTransform> & res_dirs) {
  res = st.inverse();
  double h = 1e-8;
  for (unsigned idir = 0; idir < ndirs; idir++) {
    SpatialTransform sth (
          st.E + h * st_dirs[idir].E,
          st.r + h * st_dirs[idir].r);
    SpatialTransform sth_inverse = sth.inverse();
    res_dirs[idir].E = (sth_inverse.E - res.E) / h;
    res_dirs[idir].r = (sth_inverse.r - res.r) / h;
  }
}

TEST(CheckInverse) {
  for (unsigned trial = 0; trial < 100; trial++){
    unsigned ndirs = 12u;
    Matrix3d E = Matrix3d::Random();
    Vector3d r = Vector3d::Random();
    SpatialTransform st(E, r);
    vector<SpatialTransform> st_dirs(ndirs, SpatialTransform::Zero());

    SpatialTransform ad_res = SpatialTransform::Zero();
    vector<SpatialTransform> ad_res_dirs(ndirs, SpatialTransform::Zero());

    SpatialTransform fd_res = SpatialTransform::Zero();
    vector<SpatialTransform> fd_res_dirs(ndirs, SpatialTransform::Zero());

    unsigned idir = 0;
    for (; idir < 9u; idir++) {
      st_dirs[idir].E(idir % 3, idir / 3) = 1.0;
    }
    for (; idir < 12u; idir++) {
      st_dirs[idir].r(idir - 9u) = 1.0;
    }

    AD::inverse(ndirs,
                st, st_dirs,
                ad_res, ad_res_dirs);

    FDinverse(
          ndirs,
          st, st_dirs,
          fd_res, fd_res_dirs);

    CHECK_ARRAY_CLOSE(ad_res.E.data(), fd_res.E.data(), 9, 1e-6);
    CHECK_ARRAY_CLOSE(ad_res.r.data(), fd_res.r.data(), 3, 1e-6);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      CHECK_ARRAY_CLOSE(ad_res_dirs[idir].E.data(), fd_res_dirs[idir].E.data(), 9, 1e-6);
      CHECK_ARRAY_CLOSE(ad_res_dirs[idir].r.data(), fd_res_dirs[idir].r.data(), 3, 1e-3);
    }
  }
}

// -----------------------------------------------------------------------------
