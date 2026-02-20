/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2015 Martin Felis <martin.felis@iwr.uni-heidelberg.de>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include "DynamicsAD.h"
#include "SpatialAlgebraOperatorsAD.h"

#include <iostream>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

using namespace RigidBodyDynamics::Math;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace AD {
// -----------------------------------------------------------------------------

RBDL_DLLAPI
void ForwardDynamics (Model& model,
    ADModel& ad_model,
    const VectorNd& q,
    const MatrixNd& q_dirs,
    const VectorNd& qdot,
    const MatrixNd& qdot_dirs,
    const VectorNd& tau,
    const MatrixNd& tau_dirs,
    VectorNd& qddot,
    MatrixNd& ad_qddot,
    vector<SpatialVector> const * f_ext,
    vector<vector<SpatialVector> > const * f_ext_dirs
) {
  assert((f_ext == NULL) == (f_ext_dirs == NULL));

	SpatialVector spatial_gravity (0., 0., 0., model.gravity[0], model.gravity[1], model.gravity[2]);

	unsigned int ndirs = q_dirs.cols();
	ad_model.resize_directions(ndirs);

	unsigned int i = 0;

	// Reset the velocity of the root body
	model.v[0].setZero();

	for (i = 1; i < model.mBodies.size(); i++) {
		unsigned int lambda = model.lambda[i];

		jcalc(model, ad_model, i, q, q_dirs, qdot, qdot_dirs);

		if (lambda != 0) {
//      // derivative evaluation
//      for(unsigned idir = 0; idir < ndirs; idir++) {
//        ad_model.X_base[i][idir] =
//            ad_model.X_lambda[i][idir] * model.X_base[lambda]
//            + model.X_lambda[i]/* .toMatrix()*/ * ad_model.X_base[lambda][idir];
//			}
//			// nominal evaluation
//			model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
      mulSTST(ndirs,
              model.X_lambda[i], ad_model.X_lambda[i],
              model.X_base[lambda], ad_model.X_base[lambda],
              model.X_base[i], ad_model.X_base[i]);

		} else {
			// derivative evaluation
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.X_base[i][j] = ad_model.X_lambda[i][j];
			}
			// nominal evaluation
			model.X_base[i] = model.X_lambda[i];
		}

    applySTSV(ndirs,
              model.X_lambda[i], ad_model.X_lambda[i],
              model.v[lambda], ad_model.v[lambda],
              model.v[i], ad_model.v[i]);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_model.v[i][idir] += ad_model.v_J[i][idir];
    }
    model.v[i] += model.v_J[i];

		// derivative evaluation
		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.c[i][j] = ad_model.c_J[i][j]
					+ crossm(ad_model.v[i][j], model.v_J[i])
					+ crossm(model.v[i], ad_model.v_J[i][j]);
		}
		// nominal evaluation
		model.c[i] = model.c_J[i] + crossm(model.v[i], model.v_J[i]);

		// derivative evaluation
		for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.IA[i][j].setZero();
		}
		// nominal evaluation
		model.I[i].setSpatialMatrix(model.IA[i]);

		// derivative evaluation
		for(unsigned int j = 0; j < ndirs; j++) {
			ad_model.pA[i][j] = crossf(ad_model.v[i][j],model.I[i] * model.v[i])
				+ crossf(model.v[i],model.I[i] * ad_model.v[i][j]);;
		}
		// nominal evaluation
		model.pA[i] = crossf(model.v[i],model.I[i] * model.v[i]);

    if (f_ext != NULL && (*f_ext)[i] != SpatialVector::Zero()) {
      addApplyAdjointSTSV(
            ndirs,
            -1.0,
            model.X_base[i], ad_model.X_base[i],
            (*f_ext)[i], (*f_ext_dirs)[i],
            model.pA[i], ad_model.pA[i]);
		}
	}

	for (i = model.mBodies.size() - 1; i > 0; i--) {
		unsigned int q_index = model.mJoints[i].q_index;

		if (model.mJoints[i].mDoFCount == 3) {
			cerr << "Multi-dof not supported." << endl;
			abort();
//          model.multdof3_U[i] = model.IA[i] * model.multdof3_S[i];
// #ifdef EIGEN_CORE_H
//          model.multdof3_Dinv[i] = (model.multdof3_S[i].transpose() * model.multdof3_U[i]).inverse().eval();
// #else
//          model.multdof3_Dinv[i] = (model.multdof3_S[i].transpose() * model.multdof3_U[i]).inverse();
// #endif
//          Vector3d tau_temp (Tau[q_index], Tau[q_index + 1], Tau[q_index + 2]);
//
//          model.multdof3_u[i] = tau_temp - model.multdof3_S[i].transpose() * model.pA[i];
//
// //           LOG << "multdof3_u[" << i << "] = " << model.multdof3_u[i].transpose() << std::endl;
//          unsigned int lambda = model.lambda[i];
//          if (lambda != 0) {
//              SpatialMatrix Ia = model.IA[i] - model.multdof3_U[i] * model.multdof3_Dinv[i] * model.multdof3_U[i].transpose();
//              SpatialVector pa = model.pA[i] + Ia * model.c[i] + model.multdof3_U[i] * model.multdof3_Dinv[i] * model.multdof3_u[i];
// #ifdef EIGEN_CORE_H
//              model.IA[lambda].noalias() += model.X_lambda[i].toMatrixTranspose() * Ia * model.X_lambda[i].toMatrix();
//              model.pA[lambda].noalias() += model.X_lambda[i].applyTranspose(pa);
// #else
//              model.IA[lambda] += model.X_lambda[i].toMatrixTranspose() * Ia * model.X_lambda[i].tadoMatrix();
//              model.pA[lambda] += model.X_lambda[i].applyTranspose(pa);
// #endif
//              LOG << "pA[" << lambda << "] = " << model.pA[lambda].transpose() << std::endl;
//          }
		} else {
			// derivative evaluation
			for(unsigned int j = 0; j < ndirs; j++) {
					ad_model.U[i][j] = ad_model.IA[i][j] * model.S[i];
			}
			// nominal evaluation
			model.U[i] = model.IA[i] * model.S[i];

			// derivative evaluation
			for(unsigned int j = 0; j < ndirs; j++) {
					ad_model.d(i,j) = + model.S[i].dot(ad_model.U[i][j]);
			}
			// nominal evaluation
			model.d[i] = model.S[i].dot(model.U[i]);

			// derivative evaluation
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.u(i,j) = tau_dirs(q_index,j)
						- model.S[i].dot(ad_model.pA[i][j]);
			}
			// nominal evaluation
			model.u[i] = tau[q_index] - model.S[i].dot(model.pA[i]);

			unsigned int lambda = model.lambda[i];
      if (lambda != 0) {
        vector<SpatialMatrix> ad_Ia(ndirs, SpatialMatrix::Zero());
        // derivative evaluation
				for(unsigned int j = 0; j < ndirs; j++) {
					ad_Ia[j] = ad_model.IA[i][j]
						- ad_model.U[i][j] * (model.U[i] / model.d[i]).transpose()
						- model.U[i] * (ad_model.U[i][j] / model.d[i]).transpose()
						- model.U[i] * (model.U[i] * (-ad_model.d(i,j)) / (model.d[i] * model.d[i])).transpose();
				}
				// nominal evaluation
				SpatialMatrix Ia = model.IA[i] - model.U[i] * (model.U[i] / model.d[i]).transpose();

				// derivative evaluation
				vector<SpatialVector> ad_pa(ndirs, SpatialVector::Zero());
				for(unsigned int j = 0; j < ndirs; j++) {
					ad_pa[j] = ad_model.pA[i][j]
						+ ad_Ia[j] * model.c[i]
						+ Ia * ad_model.c[i][j]
						+ ad_model.U[i][j] * model.u[i] / model.d[i]
						+ model.U[i] * ad_model.u(i,j) / model.d[i]
						+ model.U[i] * model.u[i] * (-ad_model.d(i,j)) / (model.d[i] * model.d[i]);
				}
				// nominal evaluation
				SpatialVector pa = model.pA[i] + Ia * model.c[i] + model.U[i] * model.u[i] / model.d[i];

#ifdef EIGEN_CORE_H
        addSqrFormSTSM_noalias(
              ndirs,
              model.X_lambda[i], ad_model.X_lambda[i],
              Ia, ad_Ia,
              model.IA[lambda], ad_model.IA[lambda]);

        SpatialVector summand;
        vector<SpatialVector> ad_summand(ndirs);

        applyTransposeSTSV(
              ndirs,
              model.X_lambda[i], ad_model.X_lambda[i],
              pa, ad_pa,
              summand, ad_summand);
        for (unsigned idir = 0; idir < ndirs; idir++) {
          ad_model.pA[lambda][idir].noalias() += ad_summand[idir];
        }
        model.pA[lambda].noalias() += summand;

#else
				cerr << "Simple math not yet supported." << endl;
				abort();
				//model.IA[lambda] += model.X_lambda[i].toMatrixTranspose() * Ia * model.X_lambda[i].toMatrix();
				//model.pA[lambda] += model.X_lambda[i].applyTranspose(pa);
#endif
			}
		}
	}

  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_model.a[0][idir].setZero();
  }
  model.a[0] = spatial_gravity * -1.;

	for (i = 1; i < model.mBodies.size(); i++) {
		unsigned int q_index = model.mJoints[i].q_index;
		// unsigned int lambda = model.lambda[i];
		// SpatialTransform X_lambda = model.X_lambda[i];

    SpatialVector temp = model.X_lambda[i].apply(model.a[model.lambda[i]]);

    model.a[i] = temp + model.c[i];

    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_model.a[i][idir] = crossm(temp)*model.S[i]*q_dirs(q_index, idir)
          + model.X_lambda[i].toMatrix()*ad_model.a[model.lambda[i]][idir]
          + ad_model.c[i][idir];
    }


//    applySTSV(ndirs,
//              model.X_lambda[i], ad_model.X_lambda[i],
//              model.a[lambda], ad_model.a[lambda],
//              model.a[i], ad_model.a[i]);
//    for (unsigned idir = 0; idir < ndirs; idir++) {
//      ad_model.a[i][idir] += ad_model.c[i][idir];
//    }
//    model.a[i] += model.c[i];

    if (model.mJoints[i].mDoFCount == 3
        && model.mJoints[i].mJointType != JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__
           << ": Multi-dof joints not supported." << endl;
      abort();
    } else if (model.mJoints[i].mDoFCount == 1
        && model.mJoints[i].mJointType != JointTypeCustom) {
      // derivative evaluation
      for(unsigned idir = 0; idir < ndirs; idir++) {
        ad_qddot(q_index, idir) =
            -(ad_model.d(i, idir) / (model.d[i] * model.d[i]))
            * (model.u[i] - model.U[i].dot(model.a[i]))
            + (1./model.d[i])
            * (ad_model.u(i, idir)
               - ad_model.U[i][idir].dot(model.a[i])
               - model.U[i].dot(ad_model.a[i][idir]));
      }
      // nominal evaluation
      qddot[q_index] = (1./model.d[i]) * (model.u[i] - model.U[i].dot(model.a[i]));

      // derivative evaluation
      for(unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.a[i][idir] = ad_model.a[i][idir] + model.S[i] * ad_qddot(q_index, idir);
      }
      // nominal evaluation
      model.a[i] = model.a[i] + model.S[i] * qddot[q_index];
    } else {
      cerr << __FILE__ << " " << __LINE__
           << ": Custom joints not supported." << endl;    abort();
      abort();
    }
  }
}

RBDL_DLLAPI
void InverseDynamics(
    Model& model,
    ADModel &ad_model,
    const VectorNd& q,
    const MatrixNd& q_dirs,
    const VectorNd& qdot,
    const MatrixNd& qdot_dirs,
    const VectorNd& qddot,
    const MatrixNd& qddot_dirs,
    VectorNd& tau,
    MatrixNd& ad_tau,
    vector<SpatialVector> const * f_ext,
    vector<vector<SpatialVector>> const * f_ext_dirs
) {
  assert((f_ext == NULL) == (f_ext_dirs == NULL));

  const unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs == qddot_dirs.cols());
  ad_model.resize_directions(ndirs);

  if (f_ext) {
    for (unsigned i = 0; i < f_ext_dirs->size(); i++) {
      assert(ndirs == (*f_ext_dirs)[i].size());
    }
  }

  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_model.v[0][idir].setZero();
  }
  model.v[0].setZero();

  for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_model.a[0][idir].setZero();
  }
	model.a[0].set (0., 0., 0.,
			-model.gravity[0],
			-model.gravity[1],
			-model.gravity[2]);

  std::vector<SpatialVector> h(model.mBodies.size());

	for (unsigned int i = 1; i < model.mBodies.size(); i++) {
		unsigned int lambda = model.lambda[i];
		unsigned int q_index = model.mJoints[i].q_index;

		// derivative evaluation
		jcalc (model, ad_model, i, q, q_dirs, qdot, qdot_dirs);
		// nominal evaluation
		// NOTE joints are already calculated in ad_jcalc
    // jcalc (model, i, Q, QDot);

		if (lambda != 0) {
			// derivative evaluation
//			for(unsigned int j = 0; j < ndirs; j++) {
//				ad_model.X_base[i][j] = ad_model.X_lambda[i][j] * model.X_base[lambda].toMatrix()
//					+ model.X_lambda[i].toMatrix() * ad_model.X_base[i][j];
//			}

      mulSTST(ndirs,
              model.X_lambda[i], ad_model.X_lambda[i],
              model.X_base[lambda], ad_model.X_base[lambda],
              model.X_base[i], ad_model.X_base[i]);
			// nominal evaluation
			model.X_base[i] = model.X_lambda[i] * model.X_base[lambda];
		}
    else {
			// derivative evaluation
			for(unsigned int j = 0; j < ndirs; j++) {
				ad_model.X_base[i][j] = ad_model.X_lambda[i][j];
			}
			// nominal evaluation
			model.X_base[i] = model.X_lambda[i];
		}

    applySTSV(ndirs,
              model.X_lambda[i], ad_model.X_lambda[i],
              model.v[lambda], ad_model.v[lambda],
              model.v[i], ad_model.v[i]);
    for (unsigned idir = 0; idir < ndirs; idir++) {
      ad_model.v[i][idir] += ad_model.v_J[i][idir];
    }
    model.v[i] += model.v_J[i];

    // derivative evaluation
    for(unsigned idir = 0; idir < ndirs; idir++) {
      ad_model.c[i][idir] = ad_model.c_J[i][idir]
          + crossm(ad_model.v[i][idir],model.v_J[i])
          + crossm(model.v[i],ad_model.v_J[i][idir]);
    }
    // nominal evaluation
    model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);

    if(model.mJoints[i].mJointType != JointTypeCustom) {
      if (model.mJoints[i].mDoFCount == 1) {
        applySTSV(ndirs,
                  model.X_lambda[i], ad_model.X_lambda[i],
                  model.a[lambda], ad_model.a[lambda],
                  model.a[i], ad_model.a[i]);
        for (unsigned idir = 0; idir < ndirs; idir++) {
          ad_model.a[i][idir] += ad_model.c[i][idir] + model.S[i] * qddot_dirs(q_index, idir);
        }
        model.a[i] += model.c[i] + model.S[i] * qddot(q_index);
      } else {
        cerr << __FILE__ << " " << __LINE__ << ":"
             << " Multi-DoF joint not supported." << endl;
        abort();
        // nominal evaluation
        model.a[i] = model.X_lambda[i].apply(model.a[lambda])
            + model.c[i]
            + model.multdof3_S[i] * Vector3d (qddot[q_index], qddot[q_index + 1], qddot[q_index + 2]);
      }
    } else {
      cerr << __FILE__ << " " << __LINE__ << ":"
           << " Custom joint not supported." << endl;
      abort();
      // nominal evaluation
      unsigned int k = model.mJoints[i].custom_joint_index;
      VectorNd customJointQDDot(model.mCustomJoints[k]->mDoFCount);
      for(unsigned int z = 0; z < model.mCustomJoints[k]->mDoFCount; ++z){
        customJointQDDot[z] = qddot[q_index+z];
      }
      model.a[i] =  model.X_lambda[i].apply(model.a[lambda])
        + model.c[i]
        + model.mCustomJoints[k]->S * customJointQDDot;
    }

		if (!model.mBodies[i].mIsVirtual) {
      // derivative code
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.h[i][idir] = model.I[i] * ad_model.v[i][idir];
      }
      // nominal code
      h[i] = model.I[i] * model.v[i];

      // derivative evaluation
      for(unsigned int j = 0; j < ndirs; j++) {
        ad_model.f[i][j] = model.I[i] * ad_model.a[i][j]
          + crossf(ad_model.v[i][j],model.I[i] * model.v[i])
          + crossf(model.v[i],model.I[i] * ad_model.v[i][j]);
      }
      // nominal evaluation
      model.f[i] = model.I[i] * model.a[i]
          + crossf(model.v[i],model.I[i] * model.v[i]);
		} else {
			// derivative evaluation
      for(unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.f[i][idir].setZero();
			}
			// nominal evaluation
			model.f[i].setZero();
		}
  }

  if (f_ext != NULL) {
    for (unsigned int i = 1; i < model.mBodies.size(); i++) {
      unsigned lambda = model.lambda[i];
      mulSTST(ndirs,
              model.X_lambda[i], ad_model.X_lambda[i],
              model.X_base[lambda], ad_model.X_base[lambda],
              model.X_base[i], ad_model.X_base[i]);
      addApplyAdjointSTSV(ndirs,
                          -1.0,
                          model.X_base[i], ad_model.X_base[i],
                          (*f_ext)[i], (*f_ext_dirs)[i],
                          model.f[i], ad_model.f[i]);
    }
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
    if(model.mJoints[i].mJointType != JointTypeCustom){
      if (model.mJoints[i].mDoFCount == 1) {
        // derivative evaluation
        for(unsigned idir = 0; idir < ndirs; idir++) {
          ad_tau(model.mJoints[i].q_index, idir)
              = model.S[i].dot(ad_model.f[i][idir]);
        }
        // nominal evaluation
        tau[model.mJoints[i].q_index] = model.S[i].dot(model.f[i]);
      } else if (model.mJoints[i].mDoFCount == 3) {
        cerr << __FILE__ << " " << __LINE__ << ":"
             << "Multi-DoF joint not supported." << endl;
        abort();
        // nominal
        tau.block<3,1>(model.mJoints[i].q_index, 0)
          = model.multdof3_S[i].transpose() * model.f[i];
      }
    } else if (model.mJoints[i].mJointType == JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__ << ":"
           << " Custom joint not supported." << endl;
      abort();
      // nominal
      unsigned int k = model.mJoints[i].custom_joint_index;
      tau.block(model.mJoints[i].q_index,0,
          model.mCustomJoints[k]->mDoFCount, 1)
        = model.mCustomJoints[k]->S.transpose() * model.f[i];
    }

    if (model.lambda[i] != 0) {
      SpatialVector summand;
      vector<SpatialVector> ad_summand(ndirs);
      applyTransposeSTSV(ndirs,
                         model.X_lambda[i], ad_model.X_lambda[i],
                         model.f[i], ad_model.f[i],
                         summand, ad_summand);
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.f[model.lambda[i]][idir] += ad_summand[idir];
      }
      model.f[model.lambda[i]] += summand;
    }
  }
}

RBDL_DLLAPI
void NonlinearEffects (
    Model &model,
    ADModel &ad_model,
    const VectorNd & q,
    const MatrixNd & q_dirs,
    const VectorNd & qdot,
    const MatrixNd & qdot_dirs,
    VectorNd & tau,
    MatrixNd & ad_tau
) {
  const unsigned ndirs = q_dirs.cols();
  assert(ndirs == qdot_dirs.cols());
  assert(ndirs <= ad_tau.cols());
  ad_model.resize_directions(ndirs);

  LOG << "-------- " << __func__ << " --------" << std::endl;

  SpatialVector spatial_gravity (0., 0., 0.,
      -model.gravity[0],
      -model.gravity[1],
      -model.gravity[2]);

  // Reset the velocity of the root body
  // derivative code
	for (unsigned idir = 0; idir < ndirs; idir++) {
    ad_model.v[0][idir].setZero();
    ad_model.a[0][idir].setZero();
  }
  // nominal code
  model.v[0].setZero();
  model.a[0] = spatial_gravity;

  for (unsigned i = 1; i < model.mJointUpdateOrder.size(); i++) {
    // derivative and nominal code
    jcalc(model, ad_model, model.mJointUpdateOrder[i], q, q_dirs, qdot, qdot_dirs);
  }

  std::vector<SpatialVector> h(model.mBodies.size(), SpatialVector::Zero());

  for (unsigned i = 1; i < model.mBodies.size(); i++) {
    unsigned lambda = model.lambda[i];

    if (lambda == 0) {
      // derivative code
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.v[i][idir] = ad_model.v_J[i][idir];
      }
      // nominal code
      model.v[i] = model.v_J[i];

      applySTSV(ndirs,
                model.X_lambda[i], ad_model.X_lambda[i],
                spatial_gravity,
                model.a[i], ad_model.a[i]);
    } else {
      applySTSV(ndirs,
                model.X_lambda[i], ad_model.X_lambda[i],
                model.v[lambda], ad_model.v[lambda],
                model.v[i], ad_model.v[i]);
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.v[i][idir] += ad_model.v_J[i][idir];
      }
      model.v[i] += model.v_J[i];

      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.c[i][idir] = ad_model.c_J[i][idir]
            + crossm(ad_model.v[i][idir],model.v_J[i])
            + crossm(model.v[i], ad_model.v_J[i][idir]);
      }
      model.c[i] = model.c_J[i] + crossm(model.v[i],model.v_J[i]);

      applySTSV(ndirs,
                model.X_lambda[i], ad_model.X_lambda[i],
                model.a[lambda], ad_model.a[lambda],
                model.a[i], ad_model.a[i]);
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.a[i][idir] += ad_model.c[i][idir];
      }
      model.a[i] += model.c[i];

    }

    if (!model.mBodies[i].mIsVirtual) {
      // derivative code
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.h[i][idir] = model.I[i] * ad_model.v[i][idir];
      }
      // nominal code
      h[i] = model.I[i] * model.v[i];

      // derivative code
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.f[i][idir] = model.I[i] * ad_model.a[i][idir]
            + crossf(ad_model.v[i][idir], h[i])
            + crossf(model.v[i], ad_model.h[i][idir]);
      }
      // nominal code
      model.f[i] = model.I[i] * model.a[i]
            + crossf(model.v[i], h[i]);

    } else {
      // derivative code
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.f[i][idir].setZero();
      }
      // nominal code
      model.f[i].setZero();
    }
  }

  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
    if (model.mJoints[i].mDoFCount == 3) {
      cerr << __FILE__ << " " << __LINE__ << ":"
           << "Multi-DoF joint not supported." << endl;
      abort();
      // nominal code
      tau.block<3,1>(model.mJoints[i].q_index, 0) = model.multdof3_S[i].transpose() * model.f[i];
    } else if (model.mJoints[i].mJointType == JointTypeCustom) {
      cerr << __FILE__ << " " << __LINE__ << ":"
           << " Custom joint not supported." << endl;
      abort();
    } else {
      int q_index = model.mJoints[i].q_index;
      // derivative code
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_tau(q_index, idir) = model.S[i].dot(ad_model.f[i][idir]);
      }
      // nominal code
      tau[q_index] = model.S[i].dot(model.f[i]);
    }

    unsigned lambda = model.lambda[i];
    if (lambda != 0) {
      addApplyTransposeSTSV(ndirs, 1.0,
                            model.X_lambda[i], ad_model.X_lambda[i],
                            model.f[i], ad_model.f[i],
                            model.f[lambda], ad_model.f[lambda]);
    }
  }
}

RBDL_DLLAPI void CompositeRigidBodyAlgorithm (
  Model &model,
  ADModel &ad_model,
  const VectorNd &q,
  const MatrixNd &q_dirs,
  MatrixNd &H,
  std::vector<MatrixNd> &H_ad,
  bool update_kinematics
) {

  assert (H.rows() == model.dof_count && H.cols() == model.dof_count);

  size_t ndirs = q_dirs.cols();
  ad_model.resize_directions(ndirs);

  for (unsigned int i = 1; i < model.mBodies.size(); i++) {
    if (update_kinematics) {
      // derivative evaluation
      jcalc_X_lambda_S (model, ad_model, i, q, q_dirs);
      // nominal evaluation
      // NOTE nominal evaluation is part of ad_jcalc_X_lambda_S
      // jcalc_X_lambda_S (model, i, q);
    }
    // derivative evaluation
    for (size_t idir = 0; idir < ndirs; idir++) {
      // NOTE Ic is a constant transformation wrt. q
      ad_model.Ic[i][idir].setZero();
    }
    // nominal evaluation
    model.Ic[i] = model.I[i];
  }


  for (unsigned int i = model.mBodies.size() - 1; i > 0; i--) {
    assert( model.mJoints[i].mJointType != JointTypeCustom);
    assert( model.mJoints[i].mDoFCount == 1);
    if (model.lambda[i] != 0) {
      // derivative evaluation
      // SpatialRigidBodyInertia         summand;
      // vector<SpatialRigidBodyInertia> summand_res(ndirs);
      applyTransposeSTSI(
            ndirs,
            model.X_lambda[i],
            ad_model.X_lambda[i],
            model.Ic[i],
            ad_model.Ic[i],
            ad_model.I0,
            ad_model.I_temp
            );
      for (unsigned idir = 0; idir < ndirs; idir++) {
        ad_model.Ic[model.lambda[i]][idir] += ad_model.I_temp[idir];
      }
      // nominal evaluation
      model.Ic[model.lambda[i]] += ad_model.I0;
          // + model.X_lambda[i].applyTranspose(model.Ic[i]);
    }

    unsigned int dof_index_i = model.mJoints[i].q_index;

    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; idir++) {
      // TODO this assumes 1DoF joints all over
      ad_model.F[i][idir] = ad_model.Ic[i][idir] * model.S[i];
    }
    // nominal evaluation
    SpatialVector F = model.Ic[i] * model.S[i];

    // derivative evaluation
    for (unsigned idir = 0; idir < ndirs; idir++) {
      H_ad[idir](dof_index_i, dof_index_i)
        = model.S[i].dot(ad_model.F[i][idir]);
    }
    // nominal evaluation
    H(dof_index_i, dof_index_i) = model.S[i].dot(F);

    unsigned int j = i;
    unsigned int dof_index_j = dof_index_i;


    while (model.lambda[j] != 0) {
      // // derivative evaluation
      applyTransposeSTSV(ndirs,
                         model.X_lambda[j], ad_model.X_lambda[j],
                         F, ad_model.F[i],
                         F, ad_model.F[i]);

      j = model.lambda[j];
      dof_index_j = model.mJoints[j].q_index;

      // derivative evaluation
      for (unsigned idir = 0; idir < ndirs; idir++) {
        H_ad[idir](dof_index_i,dof_index_j) = ad_model.F[i][idir].dot(model.S[j]);
        H_ad[idir](dof_index_j,dof_index_i) = H_ad[idir](dof_index_i,dof_index_j);
      }
      // nominal evaluation
      H(dof_index_i,dof_index_j) = F.dot(model.S[j]);
      H(dof_index_j,dof_index_i) = H(dof_index_i,dof_index_j);
    }

  }

}

// -----------------------------------------------------------------------------
} /* AD */
// -----------------------------------------------------------------------------
} /* RigidBodyDynamics */
// -----------------------------------------------------------------------------
