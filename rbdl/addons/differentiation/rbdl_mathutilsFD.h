/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#ifndef RBDL_MATH_FD_H
#define RBDL_MATH_FD_H

#include <cmath>
#include <limits>

#include <iostream>
#include <assert.h>

#include <rbdl/Model.h>
#include <rbdl/rbdl_mathutils.h>
#include <rbdl/SpatialAlgebraOperators.h>

#include "rbdl/Logging.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace Math {
// -----------------------------------------------------------------------------
namespace FD {
// -----------------------------------------------------------------------------

inline SpatialVector crossm (
        const SpatialVector &v1, const SpatialVector &v1_dirs,
        const SpatialVector &v2, const SpatialVector &v2_dirs
) {
    double h = 1.0e-8;

    // evaluate y(t+h*d)
    SpatialVector res = RigidBodyDynamics::Math::crossm (v1, v2);

    // evaluate y(t+h*d)
    SpatialVector res_hd = RigidBodyDynamics::Math::crossm (v1 + h * v1_dirs, v2 + h * v2_dirs);

    return (res_hd - res) / h;
}

inline SpatialMatrix crossm (
        const SpatialVector &v, const SpatialVector &v_dirs
) {
    double h = 1.0e-8;

    // evaluate y(t+h*d)
    SpatialMatrix res = RigidBodyDynamics::Math::crossm (v);

    // evaluate y(t+h*d)
    SpatialMatrix res_hd = RigidBodyDynamics::Math::crossm (v + h * v_dirs);

    return (res_hd - res) / h;
}

inline SpatialVector crossf (
        const SpatialVector &v1, const SpatialVector &v1_dirs,
        const SpatialVector &v2, const SpatialVector &v2_dirs
) {
    double h = 1.0e-8;

    // evaluate y(t+h*d)
    SpatialVector res = RigidBodyDynamics::Math::crossf (v1, v2);

    // evaluate y(t+h*d)
    SpatialVector res_hd = RigidBodyDynamics::Math::crossf (v1 + h * v1_dirs, v2 + h * v2_dirs);

    return (res_hd - res) / h;
}

inline SpatialMatrix crossf (
        const SpatialVector &v, const SpatialVector &v_dirs
) {
    double h = 1.0e-8;

    // evaluate y(t+h*d)
    SpatialMatrix res = RigidBodyDynamics::Math::crossf (v);

    // evaluate y(t+h*d)
    SpatialMatrix res_hd = RigidBodyDynamics::Math::crossf (v + h * v_dirs);

    return (res_hd - res) / h;
}

inline void addApplyAdjointSTSV (
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

  SpatialVector upd   = st.toMatrixAdjoint() * sv;
  double h = 1.0e-8;
  for (unsigned idir = 0; idir < ndirs; idir++) {
    SpatialVector upd_h = SpatialTransform(st.E + h * st_dirs[idir].E, st.r + h*st_dirs[idir].r).toMatrixAdjoint()
        * (sv + h * sv_dirs[idir]);
    res_dirs[idir] += dscal * (upd_h - upd) / h;
  }
  res += dscal * upd;
}

inline void addApplyAdjointSTSV (
    unsigned ndirs,
    double dscal,
    SpatialTransform const & st,
    std::vector<SpatialTransform> const & st_dirs,
    SpatialVector const & sv,
    MatrixNd const & sv_dirs,
    SpatialVector & res,
    std::vector<SpatialVector> & res_dirs) {
  assert(ndirs <= st_dirs.size());
  assert(ndirs <= sv_dirs.cols());
  assert(ndirs <= res_dirs.size());

  SpatialVector upd   = st.toMatrixAdjoint() * sv;
  double h = 1.0e-8;
  for (unsigned idir = 0; idir < ndirs; idir++) {
    SpatialVector upd_h = SpatialTransform(st.E + h * st_dirs[idir].E, st.r + h*st_dirs[idir].r).toMatrixAdjoint()
        * (sv + h * sv_dirs.col(idir));
    res_dirs[idir] += dscal * (upd_h - upd) / h;
  }
  res += dscal * upd;
}

// -----------------------------------------------------------------------------
} /* FD */
// -----------------------------------------------------------------------------
} /* Math */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------

/* RBDL_MATH_FD_H */
#endif
