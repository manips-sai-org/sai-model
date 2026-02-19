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

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
Vector3d fd_jcalc(
    Model &model,
    unsigned int joint_id,
    const VectorNd &q,
    const MatrixNd &q_dirs,
    const VectorNd &qdot,
    const MatrixNd &qdot_dirs,
    std::vector<SpatialMatrix> &ad_X_Ji,
    std::vector<SpatialVector> &ad_S_i,
    std::vector<SpatialVector> &ad_v_Ji,
    std::vector<SpatialVector> &ad_c_Ji
) {
    double h = 1.0e-8;
    unsigned int ndirs = q_dirs.cols();

    // check input dimensions
    if (q_dirs.cols() != qdot_dirs.cols()) {
        std::cerr << "directions have different dimensions: " << "#q_dirs = " << q_dirs.cols() << " != " << qdot_dirs.cols() << " = #qdot_dirs." << std::endl;
        std::cerr << "In: " << __func__ << endl;
        abort();
    }
    // check output dimensions
    if (ad_X_Ji.size() != ndirs) {
        std::cerr << "derivative does not have proper dimensions " << "#ad_X_Ji.size() = " << ad_X_Ji.size() << " != " << ndirs << " = ndirs." << std::endl;
        std::cerr << "In: " << __func__ << endl;
        abort();
    }
    if (ad_S_i.size() != ndirs) {
        std::cerr << "derivative does not have proper dimensions " << "#ad_S_i.size() = " << ad_S_i.size() << " != " << ndirs << " = ndirs." << std::endl;
        std::cerr << "In: " << __func__ << endl;
        abort();
    }
    if (ad_c_Ji.size() != ndirs) {
        std::cerr << "derivative does not have proper dimensions " << "#ad_c_Ji.size() = " << ad_c_Ji.size() << " != " << ndirs << " = ndirs." << std::endl;
        std::cerr << "In: " << __func__ << endl;
        abort();
    }

    // calculate y(t)
    jcalc (model, joint_id, q, qdot);
    MatrixNd ref_X_Ji = model.X_J[joint_id].toMatrix();
    SpatialVector ref_S_i = model.S[joint_id];
    SpatialVector ref_v_Ji = model.v_J[joint_id];
    SpatialVector ref_c_Ji = model.c_J[joint_id];

    for (unsigned int j = 0; j < ndirs; j++) {
        VectorNd q_dir = q_dirs.block(0,j, model.qdot_size, 1);
        VectorNd qdot_dir = qdot_dirs.block(0,j, model.qdot_size, 1);
        jcalc (model, joint_id, q + h * q_dir, qdot + h * qdot_dir);

        MatrixNd hd_X_Ji = model.X_J[joint_id].toMatrix();
        SpatialVector hd_S_i = model.S[joint_id];
        SpatialVector hd_v_Ji = model.v_J[joint_id];
        SpatialVector hd_c_Ji = model.c_J[joint_id];

        ad_X_Ji[j] = (hd_X_Ji - ref_X_Ji) / h;
        ad_S_i[j]  = (hd_S_i  - ref_S_i)  / h;
        ad_v_Ji[j] = (hd_v_Ji - ref_v_Ji) / h;
        ad_c_Ji[j] = (hd_c_Ji - ref_c_Ji) / h;
    }
};

// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------
