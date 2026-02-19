#include <cfloat>
#include <cmath>

#include "rbdl_utilsFD.h"

#include <rbdl/rbdl_utils.h>
#include <rbdl/Model.h>

#include "FdModelEntry.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace Utils {
// -----------------------------------------------------------------------------
namespace FDC {
// -----------------------------------------------------------------------------

using namespace Math;

RBDL_DLLAPI void CalcCenterOfMass (
    Model & model,
    ADModel * fd_model,
    VectorNd const & q,
    MatrixNd const & q_dirs,
    VectorNd const & qdot,
    MatrixNd const & qdot_dirs,
    double & mass,
    Vector3d & com,
    MatrixNd & fd_com,
    Vector3d * com_velocity,
    MatrixNd * fd_com_velocity,
    Vector3d * angular_momentum,
    MatrixNd * fd_angular_momentum
) {
  size_t ndirs = q_dirs.cols();
  assert(ndirs == static_cast<unsigned>(qdot_dirs.cols()));
  assert(ndirs == static_cast<unsigned>(fd_com.cols()));

  if (com_velocity) {
    assert(q_dirs.cols() == fd_com_velocity->cols());
  }
  if (angular_momentum) {
    assert(q_dirs.cols() == fd_angular_momentum->cols());
  }

  RigidBodyDynamics::Utils::CalcCenterOfMass(model, q, qdot, mass, com,
      com_velocity, angular_momentum, true);

  VectorNd q_dir(q);
  VectorNd qdot_dir(qdot);

  for (size_t i = 0; i < ndirs; i++ ) {
    Model * modelh;
    if (fd_model) {
      modelh = new Model(model);
    } else {
      modelh = &model;
    }
    q_dir = q_dirs.block(0, i, model.q_size, 1);
    qdot_dir = qdot_dirs.block(0, i, model.qdot_size, 1);

    double hd_mass;
    Vector3d hd_com;
    Vector3d * hd_com_velocity = 0;
    Vector3d * hd_angular_momentum = 0;
    if (com_velocity) {
      hd_com_velocity = new Vector3d;
    }
    if (angular_momentum) {
      hd_angular_momentum = new Vector3d;
    }

    Utils::CalcCenterOfMass(*modelh, q + h * q_dir,
      qdot + h * qdot_dir, hd_mass, hd_com,
      hd_com_velocity, hd_angular_momentum, true);

    fd_com.block<3,1>(0, i) = (hd_com - com) / h;
    if (com_velocity) {
      fd_com_velocity->block<3,1>(0, i) =
          (*hd_com_velocity - *com_velocity) / h;
      delete hd_com_velocity;
    }
    if (angular_momentum) {
      fd_angular_momentum->block<3,1>(0, i) =
          (*hd_angular_momentum - *angular_momentum) / h;
      delete hd_angular_momentum;
    }

    if (fd_model) {
      computeFDEntry(model, *modelh, h, i, *fd_model);
      delete modelh;
    }
  }
}

RBDL_DLLAPI double CalcPotentialEnergy (
        Model & model,
        VectorNd const & q,
        MatrixNd const & q_dirs,
        MatrixNd & fd_pote) {
    unsigned int ndirs = q_dirs.cols();

    assert(ndirs == fd_pote.cols());
    assert(1 == fd_pote.rows());

    double pote = RigidBodyDynamics::Utils::CalcPotentialEnergy(model, q);

    VectorNd q_dir(VectorNd::Zero(ndirs));
    for (unsigned int i = 0; i < ndirs; i++ ) {
        q_dir = q_dirs.block(0, i, model.q_size, 1);
        double hd_pote = RigidBodyDynamics::Utils::CalcPotentialEnergy(model,
                q + h * q_dir);
        fd_pote(i) = (hd_pote - pote) / h;
    }

    return pote;
}

RBDL_DLLAPI double CalcKineticEnergy (
        Model & model,
        VectorNd const & q,
        MatrixNd const & q_dirs,
        VectorNd const & qdot,
        MatrixNd const & qdot_dirs,
        MatrixNd & fd_kine) {
    unsigned int ndirs = q_dirs.cols();

    assert(ndirs == qdot_dirs.cols());
    assert(ndirs == fd_kine.cols());
    assert(1 == fd_kine.rows());

    double kine = RigidBodyDynamics::Utils::CalcKineticEnergy(model, q, qdot);

    VectorNd q_dir(VectorNd::Zero(ndirs));
    VectorNd qdot_dir(VectorNd::Zero(ndirs));
    for (unsigned int i = 0; i < ndirs; i++ ) {
        q_dir = q_dirs.block(0, i, model.q_size, 1);
        qdot_dir = qdot_dirs.block(0, i, model.qdot_size, 1);
        double hd_kine = RigidBodyDynamics::Utils::CalcKineticEnergy(model,
                q + h * q_dir, qdot + h * qdot_dir);
        fd_kine(i) = (hd_kine - kine) / h;
    }

    return kine;
}

// -----------------------------------------------------------------------------
} // namespace FD
// -----------------------------------------------------------------------------
} // namespace Utils
// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------
