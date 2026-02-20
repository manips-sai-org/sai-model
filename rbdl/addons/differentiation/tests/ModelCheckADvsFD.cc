#include <UnitTest++.h>

#include "ModelCheckADvsFD.h"

using namespace std;

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

void checkModelsADvsFD(
    unsigned ndirs,
    Model const & ad_model,
    ADModel const & ad_d_model,
    Model const & fd_model,
    ADModel const & fd_d_model,
    double const& PREC
) {

  // nominal check
  CHECK_CLOSE(ad_model.v_J.size(), fd_model.v_J.size(), PREC);
  for (unsigned i = 0; i < ad_model.v_J.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.v_J[i].data(), fd_model.v_J[i].data(), 6, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.v_J.size(); i++) {
      if((ad_d_model.v_J[i][idir] - fd_d_model.v_J[i][idir]).norm() > PREC) {
        std::cout << "v_J " << i << "," << idir << std::endl;
        std::cout << ad_d_model.v_J[i][idir].transpose() << std::endl;
        std::cout << fd_d_model.v_J[i][idir].transpose() << std::endl;
      }
      CHECK_ARRAY_CLOSE(ad_d_model.v_J[i][idir].data(), fd_d_model.v_J[i][idir].data(), 6, PREC );
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.c_J.size(), fd_model.c_J.size(), PREC);
  for (unsigned i = 0; i < ad_model.c_J.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.c_J[i].data(), fd_model.c_J[i].data(), 6, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.c_J.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.c_J[i][idir].data(), fd_d_model.c_J[i][idir].data(), 6, PREC );
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.X_J.size(), fd_model.X_J.size(), PREC);
  for (unsigned i = 0; i < ad_model.X_J.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.X_J[i].E.data(), fd_model.X_J[i].E.data(), 9, PREC);
    CHECK_ARRAY_CLOSE(ad_model.X_J[i].r.data(), fd_model.X_J[i].r.data(), 3, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.X_J.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.X_J[i][idir].E.data(), fd_d_model.X_J[i][idir].E.data(), 9, PREC);
      CHECK_ARRAY_CLOSE(ad_d_model.X_J[i][idir].r.data(), fd_d_model.X_J[i][idir].r.data(), 3, PREC);
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.Ic.size(), fd_model.Ic.size(), PREC);
  for (unsigned i = 0; i < ad_model.Ic.size(); i++) {
    CHECK_CLOSE(ad_model.Ic[i].m, fd_model.Ic[i].m, PREC);
    CHECK_ARRAY_CLOSE(ad_model.Ic[i].h.data(), fd_model.Ic[i].h.data(), 3, PREC);
    CHECK_CLOSE(ad_model.Ic[i].Ixx, fd_model.Ic[i].Ixx, PREC);
    CHECK_CLOSE(ad_model.Ic[i].Iyx, fd_model.Ic[i].Iyx, PREC);
    CHECK_CLOSE(ad_model.Ic[i].Iyy, fd_model.Ic[i].Iyy, PREC);
    CHECK_CLOSE(ad_model.Ic[i].Izx, fd_model.Ic[i].Izx, PREC);
    CHECK_CLOSE(ad_model.Ic[i].Izy, fd_model.Ic[i].Izy, PREC);
    CHECK_CLOSE(ad_model.Ic[i].Izz, fd_model.Ic[i].Izz, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.Ic.size(); i++) {
      CHECK_CLOSE(ad_d_model.Ic[i][idir].m, fd_d_model.Ic[i][idir].m, PREC);
      CHECK_ARRAY_CLOSE(ad_d_model.Ic[i][idir].h.data(), fd_d_model.Ic[i][idir].h.data(), 3, PREC);
      CHECK_CLOSE(ad_d_model.Ic[i][idir].Ixx, fd_d_model.Ic[i][idir].Ixx, PREC);
      CHECK_CLOSE(ad_d_model.Ic[i][idir].Iyx, fd_d_model.Ic[i][idir].Iyx, PREC);
      CHECK_CLOSE(ad_d_model.Ic[i][idir].Iyy, fd_d_model.Ic[i][idir].Iyy, PREC);
      CHECK_CLOSE(ad_d_model.Ic[i][idir].Izx, fd_d_model.Ic[i][idir].Izx, PREC);
      CHECK_CLOSE(ad_d_model.Ic[i][idir].Izy, fd_d_model.Ic[i][idir].Izy, PREC);
      CHECK_CLOSE(ad_d_model.Ic[i][idir].Izz, fd_d_model.Ic[i][idir].Izz, PREC);
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.X_lambda.size(), fd_model.X_lambda.size(), PREC);
  for (unsigned i = 0; i < ad_model.X_lambda.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.X_lambda[i].E.data(), fd_model.X_lambda[i].E.data(), 9, PREC);
    CHECK_ARRAY_CLOSE(ad_model.X_lambda[i].r.data(), fd_model.X_lambda[i].r.data(), 3, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.X_lambda.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.X_lambda[i][idir].E.data(), fd_d_model.X_lambda[i][idir].E.data(), 9, PREC);
      CHECK_ARRAY_CLOSE(ad_d_model.X_lambda[i][idir].r.data(), fd_d_model.X_lambda[i][idir].r.data(), 3, PREC);
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.U.size(), fd_model.U.size(), PREC);
  for (unsigned i = 0; i < ad_model.U.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.U[i].data(), fd_model.U[i].data(), 6, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.U.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.U[i][idir].data(), fd_d_model.U[i][idir].data(), 6, PREC);
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.d.rows(), fd_model.d.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_model.d.data(), fd_model.d.data(), ad_model.d.rows(), PREC);
  // derivative check
  CHECK_ARRAY_CLOSE(ad_d_model.d.data(), fd_d_model.d.data(), ad_d_model.u.cols() * ad_model.d.rows(), PREC);

  // nominal check
  CHECK_CLOSE(ad_model.IA.size(), fd_model.IA.size(), PREC);
  for (unsigned i = 0; i < ad_model.IA.size(); i++) {
//    std::cout << "--- AD --- " << i << std::endl;
//    std::cout << ad_model.IA[i] << std::endl;
//    std::cout << "--- FD --- " << i << std::endl;
//    std::cout << fd_model.IA[i] << std::endl;
    CHECK_ARRAY_CLOSE(ad_model.IA[i].data(), fd_model.IA[i].data(), 36, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.IA.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.IA[i][idir].data(), fd_d_model.IA[i][idir].data(), 36, PREC);
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.X_base.size(), fd_model.X_base.size(), PREC);
  for (unsigned i = 0; i < ad_model.X_base.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.X_base[i].E.data(), fd_model.X_base[i].E.data(), 9, PREC);
    CHECK_ARRAY_CLOSE(ad_model.X_base[i].r.data(), fd_model.X_base[i].r.data(), 3, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.X_base.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.X_base[i][idir].E.data(), fd_d_model.X_base[i][idir].E.data(), 9, PREC);
      CHECK_ARRAY_CLOSE(ad_d_model.X_base[i][idir].r.data(), fd_d_model.X_base[i][idir].r.data(), 3, PREC);
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.v.size(), fd_model.v.size(), PREC);
  for (unsigned i = 0; i < ad_model.v.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.v[i].data(), fd_model.v[i].data(), 6, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.v.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.v[i][idir].data(), fd_d_model.v[i][idir].data(), 6, PREC );
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.c.size(), fd_model.c.size(), PREC);
  for (unsigned i = 0; i < ad_model.c.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.c[i].data(), fd_model.c[i].data(), 6, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.c.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.c[i][idir].data(), fd_d_model.c[i][idir].data(), 6, PREC );
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.a.size(), fd_model.a.size(), PREC);
  for (unsigned i = 0; i < ad_model.a.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.a[i].data(), fd_model.a[i].data(), 6, PREC);
  }

  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.a.size(); i++) {
      if ( (ad_d_model.a[i][idir] - fd_d_model.a[i][idir]).norm() > 5e-5) {
        std::cout << "der a " << i << "," << idir << std::endl;
        std::cout << ad_d_model.a[i][idir].transpose() << std::endl;
        std::cout << fd_d_model.a[i][idir].transpose() << std::endl;
      }
      CHECK_ARRAY_CLOSE(ad_d_model.a[i][idir].data(), fd_d_model.a[i][idir].data(), 6, 5e-5);
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.pA.size(), fd_model.pA.size(), PREC);
  for (unsigned i = 0; i < ad_model.pA.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.pA[i].data(), fd_model.pA[i].data(), 6, PREC );
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.pA.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.pA[i][idir].data(), fd_d_model.pA[i][idir].data(), 6, 5e-5);
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.u.rows(), fd_model.u.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_model.u.data(), fd_model.u.data(), ad_model.u.rows(), PREC);
  // derivative check
  CHECK_ARRAY_CLOSE(ad_d_model.u.data(), fd_d_model.u.data(), ad_d_model.u.cols() * ad_model.u.rows(), 5e-5);

  // nominal check
  CHECK_CLOSE(ad_model.f.size(), fd_model.f.size(), PREC);
  for (unsigned i = 0; i < ad_model.f.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.f[i].data(), fd_model.f[i].data(), 6, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.f.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.f[i][idir].data(), fd_d_model.f[i][idir].data(), 6, 5e-5);
    }
  }
}


void checkModelsADvsED(
    unsigned ndirs,
    Model const & ad_model,
    ADModel const & ad_d_model,
    Model const & ed_model,
    EDModel const & ed_d_model,
    double const& PREC
) {

  // nominal check
  CHECK_CLOSE(ad_model.v_J.size(), ed_model.v_J.size(), PREC);
  for (unsigned i = 0; i < ad_model.v_J.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.v_J[i].data(), ed_model.v_J[i].data(), 6, PREC);
  }
  // // derivative check
//   for (unsigned idir = 0; idir < ndirs; idir++) {
//     for (unsigned i = 0; i < ad_model.v_J.size(); i++) {
//       CHECK_ARRAY_CLOSE(ad_d_model.v_J[i][idir].data(), ed_d_model.v_J[i].col(idir).data(), 6, PREC );
//     }
//   }

  // nominal check
  CHECK_CLOSE(ad_model.c_J.size(), ed_model.c_J.size(), PREC);
  for (unsigned i = 0; i < ad_model.c_J.size(); i++) {   
    CHECK_ARRAY_CLOSE(ad_model.c_J[i].data(), ed_model.c_J[i].data(), 6, PREC);
  }
  // derivative check
//   for (unsigned idir = 0; idir < ndirs; idir++) {
//     for (unsigned i = 0; i < ad_model.c_J.size(); i++) {
//       CHECK_ARRAY_CLOSE(ad_d_model.c_J[i][idir].data(), ed_d_model.c_J[i].col(idir).data(), 6, PREC );
//     }
//   }

  // nominal check
  CHECK_CLOSE(ad_model.X_J.size(), ed_model.X_J.size(), PREC);
  for (unsigned i = 0; i < ad_model.X_J.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.X_J[i].E.data(), ed_model.X_J[i].E.data(), 9, PREC);
    CHECK_ARRAY_CLOSE(ad_model.X_J[i].r.data(), ed_model.X_J[i].r.data(), 3, PREC);
  }
  // derivative check
  // for (unsigned idir = 0; idir < ndirs; idir++) {
  //   for (unsigned i = 0; i < ad_model.X_J.size(); i++) {
  //     CHECK_ARRAY_CLOSE(ad_d_model.X_J[i][idir].E.data(), ed_d_model.X_J[i][idir].E.data(), 9, PREC);
  //     CHECK_ARRAY_CLOSE(ad_d_model.X_J[i][idir].r.data(), ed_d_model.X_J[i][idir].r.data(), 3, PREC);
  //   }
  // }

  // derivative check
  if (true)
  {
    CHECK_CLOSE(ad_d_model.F.size(), ed_d_model.F.size(), PREC);
    for (unsigned i = 0; i < ad_d_model.F.size(); i++) {
      for (unsigned idir = 0; idir < ndirs; idir++) {
        // DEBUG OUTPUT
        RigidBodyDynamics::Math::SpatialVector error
          = (ad_d_model.F[i][idir] - ed_d_model.F[i].col(idir)).cwiseAbs();
        double max = error.maxCoeff();
        if (max > PREC) {
          std::cout << "error [" << i << "][" << idir << "] = \n" << error << std::endl;
          std::cout << "max   [" << i << "][" << idir << "] = " << max << std::endl;
          std::cout << "ad_d_model.F[" << i << "][" << idir << "] = \n"
            << ad_d_model.F[i][idir] << std::endl;
          std::cout << "ed_d_model.F[" << i << "][" << idir << "] = \n"
            << ed_d_model.F[i].col(idir) << std::endl;
          std::cout << endl;
        }
        CHECK_ARRAY_CLOSE(
          ad_d_model.F[i][idir].data(),
          ed_d_model.F[i].col(idir).data(),
          6, PREC
        );
      }
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.Ic.size(), ed_model.Ic.size(), PREC);
  for (unsigned i = 0; i < ad_model.Ic.size(); i++) {
    CHECK_CLOSE(ad_model.Ic[i].m, ed_model.Ic[i].m, PREC);
    CHECK_ARRAY_CLOSE(ad_model.Ic[i].h.data(), ed_model.Ic[i].h.data(), 3, PREC);
    CHECK_CLOSE(ad_model.Ic[i].Ixx, ed_model.Ic[i].Ixx, PREC);
    CHECK_CLOSE(ad_model.Ic[i].Iyx, ed_model.Ic[i].Iyx, PREC);
    CHECK_CLOSE(ad_model.Ic[i].Iyy, ed_model.Ic[i].Iyy, PREC);
    CHECK_CLOSE(ad_model.Ic[i].Izx, ed_model.Ic[i].Izx, PREC);
    CHECK_CLOSE(ad_model.Ic[i].Izy, ed_model.Ic[i].Izy, PREC);
    CHECK_CLOSE(ad_model.Ic[i].Izz, ed_model.Ic[i].Izz, PREC);
  }
  // derivative check
  for (unsigned i = 0; i < ad_model.Ic.size(); i++) {
    for (unsigned idir = 0; idir < ndirs; idir++) {
      // DEBUG OUTPUT
      RigidBodyDynamics::Math::SpatialMatrix error
        = (ad_d_model.Ic[i][idir].toMatrix() - ed_d_model.Ic[i][idir].toMatrix()).cwiseAbs();
      double max = error.maxCoeff();
      if (max > PREC) {
        std::cout << "error [" << i << "][" << idir << "] = \n" << error << std::endl;
        std::cout << "max   [" << i << "][" << idir << "] = " << max << std::endl;
        std::cout << "ad_d_model.Ic[" << i << "][" << idir << "] = \n"
          << ad_d_model.Ic[i][idir].toMatrix() << std::endl;
        std::cout << "ed_d_model.Ic[" << i << "][" << idir << "] = \n"
          << ed_d_model.Ic[i][idir].toMatrix() << std::endl;
        std::cout << endl;
      }
      CHECK_ARRAY_CLOSE(
        ad_d_model.Ic[i][idir].toMatrix().data(), ed_d_model.Ic[i][idir].toMatrix().data(),
        6*6, PREC
      );
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.X_lambda.size(), ed_model.X_lambda.size(), PREC);
  for (unsigned i = 0; i < ad_model.X_lambda.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.X_lambda[i].E.data(), ed_model.X_lambda[i].E.data(), 9, PREC);
    CHECK_ARRAY_CLOSE(ad_model.X_lambda[i].r.data(), ed_model.X_lambda[i].r.data(), 3, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.X_lambda.size(); i++) {
      // DEBUG OUTPUT
      RigidBodyDynamics::Math::Matrix3d X_E_error
        = (ad_d_model.X_lambda[i][idir].E - ed_d_model.X_lambda[i][idir].E).cwiseAbs();
      RigidBodyDynamics::Math::Vector3d X_r_error
        = (ad_d_model.X_lambda[i][idir].r - ed_d_model.X_lambda[i][idir].r).cwiseAbs();
      const double error = X_E_error.maxCoeff() + X_r_error.maxCoeff();
      if (error > PREC) {
        std::cout << "error   [" << idir << "] = " << error << std::endl;
        std::cout << "E_error [" << idir << "] = \n" << X_E_error << std::endl;
        std::cout << "ad_d_model.X_lambda[" << i << "][" << idir << "].E = \n"
          << ad_d_model.X_lambda[i][idir].E << std::endl;
        std::cout << "ed_d_model.X_lambda[" << i << "][" << idir << "].E = \n"
          << ed_d_model.X_lambda[i][idir].E << std::endl;
        std::cout << "r_error [" << idir << "] = \n" << X_r_error.transpose() << std::endl;
        std::cout << "ad_d_model.X_lambda[" << i << "][" << idir << "].r = \n"
          << ad_d_model.X_lambda[i][idir].r << std::endl;
        std::cout << "ed_d_model.X_lambda[" << i << "][" << idir << "].r = \n"
          << ed_d_model.X_lambda[i][idir].r << std::endl;
        std::cout << endl;
      }

      CHECK_ARRAY_CLOSE(ad_d_model.X_lambda[i][idir].E.data(), ed_d_model.X_lambda[i][idir].E.data(), 9, PREC);
      CHECK_ARRAY_CLOSE(ad_d_model.X_lambda[i][idir].r.data(), ed_d_model.X_lambda[i][idir].r.data(), 3, PREC);
    }
  }

  // nominal check
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
      // DEBUG OUTPUT
      RigidBodyDynamics::Math::Matrix3d X_E_error
        = (ad_d_model.X_0[idir].E - ed_d_model.X_0[idir].E).cwiseAbs();
      RigidBodyDynamics::Math::Vector3d X_r_error
        = (ad_d_model.X_0[idir].r - ed_d_model.X_0[idir].r).cwiseAbs();
      const double error = X_E_error.maxCoeff() + X_r_error.maxCoeff();
      if (error > PREC) {
        std::cout << "error   [" << idir << "] = " << error << std::endl;
        std::cout << "E_error [" << idir << "] = \n" << X_E_error << std::endl;
        std::cout << "ad_d_model.X_0[" << idir << "].E = \n"
          << ad_d_model.X_0[idir].E << std::endl;
        std::cout << "ed_d_model.X_0[" << idir << "].E = \n"
          << ed_d_model.X_0[idir].E << std::endl;
        std::cout << "r_error [" << idir << "] = \n" << X_r_error.transpose() << std::endl;
        std::cout << "ad_d_model.X_0[" << idir << "].r = \n"
          << ad_d_model.X_0[idir].r << std::endl;
        std::cout << "ed_d_model.X_0[" << idir << "].r = \n"
          << ed_d_model.X_0[idir].r << std::endl;
        std::cout << endl;
      }
    CHECK_ARRAY_CLOSE(ad_d_model.X_0[idir].E.data(), ed_d_model.X_0[idir].E.data(), 9, PREC);
    CHECK_ARRAY_CLOSE(ad_d_model.X_0[idir].r.data(), ed_d_model.X_0[idir].r.data(), 3, PREC);
  }


  /*
  // nominal check
  CHECK_CLOSE(ad_model.U.size(), ed_model.U.size(), PREC);
  for (unsigned i = 0; i < ad_model.U.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.U[i].data(), ed_model.U[i].data(), 6, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.U.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.U[i][idir].data(), ed_d_model.U[i].col(idir).data(), 6, PREC);
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.d.rows(), ed_model.d.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_model.d.data(), ed_model.d.data(), ad_model.d.rows(), PREC);
  // derivative check
  CHECK_ARRAY_CLOSE(ad_d_model.d.data(), ed_d_model.d.data(), ad_d_model.u.cols() * ad_model.d.rows(), PREC);

  // nominal check
  CHECK_CLOSE(ad_model.IA.size(), ed_model.IA.size(), PREC);
  for (unsigned i = 0; i < ad_model.IA.size(); i++) {
//    std::cout << "--- AD --- " << i << std::endl;
//    std::cout << ad_model.IA[i] << std::endl;
//    std::cout << "--- FD --- " << i << std::endl;
//    std::cout << ed_model.IA[i] << std::endl;
    CHECK_ARRAY_CLOSE(ad_model.IA[i].data(), ed_model.IA[i].data(), 36, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.IA.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.IA[i][idir].data(), ed_d_model.IA[i][idir].data(), 36, PREC);
    }
  }
*/

  // nominal check
  // std::cout << "Check X_base disabled!" << std::endl;
  CHECK_CLOSE(ad_model.X_base.size(), ed_model.X_base.size(), PREC);
  for (unsigned i = 0; i < ad_model.X_base.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.X_base[i].E.data(), ed_model.X_base[i].E.data(), 9, PREC);
    CHECK_ARRAY_CLOSE(ad_model.X_base[i].r.data(), ed_model.X_base[i].r.data(), 3, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.X_base.size(); i++) {
      RigidBodyDynamics::Math::Matrix3d X_E_error
        = (ad_d_model.X_base[i][idir].E - ed_d_model.X_base[i][idir].E).cwiseAbs();
      RigidBodyDynamics::Math::Vector3d X_r_error
        = (ad_d_model.X_base[i][idir].r - ed_d_model.X_base[i][idir].r).cwiseAbs();
      const double error = X_E_error.maxCoeff() + X_r_error.maxCoeff();
      if (error > PREC) {
        std::cout << "error   [" << idir << "] = " << error << std::endl;
        std::cout << "E_error [" << idir << "] = \n" << X_E_error << std::endl;
        std::cout << "ad_d_model.X_base[" << i << "][" << idir << "].E = \n"
          << ad_d_model.X_base[i][idir].E << std::endl;
        std::cout << "ed_d_model.X_base[" << i << "][" << idir << "].E = \n"
          << ed_d_model.X_base[i][idir].E << std::endl;
        std::cout << "r_error [" << idir << "] = \n" << X_r_error.transpose() << std::endl;
        std::cout << "ad_d_model.X_base[" << i << "][" << idir << "].r = \n"
          << ad_d_model.X_base[i][idir].r << std::endl;
        std::cout << "ed_d_model.X_base[" << i << "][" << idir << "].r = \n"
          << ed_d_model.X_base[i][idir].r << std::endl;
        std::cout << endl;
      }

      CHECK_ARRAY_CLOSE(ad_d_model.X_base[i][idir].E.data(), ed_d_model.X_base[i][idir].E.data(), 9, PREC);
      CHECK_ARRAY_CLOSE(ad_d_model.X_base[i][idir].r.data(), ed_d_model.X_base[i][idir].r.data(), 3, PREC);
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.v.size(), ed_model.v.size(), PREC);
  for (unsigned i = 0; i < ad_model.v.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.v[i].data(), ed_model.v[i].data(), 6, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.v.size(); i++) {
      RigidBodyDynamics::Math::SpatialVector v_error
        = (ad_d_model.v[i][idir] - ed_d_model.v[i].col(idir)).cwiseAbs();
      const double error = v_error.maxCoeff();
      if (error > PREC) {
        std::cout << "error       [" << i << "][" << idir << "] = " << error << std::endl;
        std::cout << "v_error     [" << i << "][" << idir << "] = " << v_error.transpose() << std::endl;
        std::cout << "ad_d_model.v[" << i << "][" << idir << "] = "
          << ad_d_model.v[i][idir].transpose() << std::endl;
        std::cout << "ed_d_model.v[" << i << "][" << idir << "] = "
          << ed_d_model.v[i].col(idir).transpose() << std::endl;
        std::cout << endl;
      }
      CHECK_ARRAY_CLOSE(ad_d_model.v[i][idir].data(), ed_d_model.v[i].col(idir).data(), 6, PREC );
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.c.size(), ed_model.c.size(), PREC);
  for (unsigned i = 0; i < ad_model.c.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.c[i].data(), ed_model.c[i].data(), 6, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.c.size(); i++) {
      RigidBodyDynamics::Math::SpatialVector v_error
        = (ad_d_model.c[i][idir] - ed_d_model.c[i].col(idir)).cwiseAbs();
      const double error = v_error.maxCoeff();
      if (error > PREC) {
        std::cout << "error       [" << i << "][" << idir << "] = " << error << std::endl;
        std::cout << "v_error     [" << i << "][" << idir << "] = " << v_error.transpose() << std::endl;
        std::cout << "ad_d_model.c[" << i << "][" << idir << "] = "
          << ad_d_model.c[i][idir].transpose() << std::endl;
        std::cout << "ed_d_model.c[" << i << "][" << idir << "] = "
          << ed_d_model.c[i].col(idir).transpose() << std::endl;
        std::cout << endl;
      }
      CHECK_ARRAY_CLOSE(ad_d_model.c[i][idir].data(), ed_d_model.c[i].col(idir).data(), 6, PREC );
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.a.size(), ed_model.a.size(), PREC);
  for (unsigned i = 0; i < ad_model.a.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.a[i].data(), ed_model.a[i].data(), 6, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.a.size(); i++) {
      RigidBodyDynamics::Math::SpatialVector v_error
        = (ad_d_model.a[i][idir] - ed_d_model.a[i].col(idir)).cwiseAbs();
      const double error = v_error.maxCoeff();
      if (error > PREC) {
        std::cout << "error       [" << i << "][" << idir << "] = " << error << std::endl;
        std::cout << "v_error     [" << i << "][" << idir << "] = " << v_error.transpose() << std::endl;
        std::cout << "ad_d_model.a[" << i << "][" << idir << "] = "
          << ad_d_model.a[i][idir].transpose() << std::endl;
        std::cout << "ed_d_model.a[" << i << "][" << idir << "] = "
          << ed_d_model.a[i].col(idir).transpose() << std::endl;
        std::cout << endl;
      }
      CHECK_ARRAY_CLOSE(ad_d_model.a[i][idir].data(), ed_d_model.a[i].col(idir).data(), 6, 5e-5);
    }
  }

/*
  // nominal check
  CHECK_CLOSE(ad_model.pA.size(), ed_model.pA.size(), PREC);
  for (unsigned i = 0; i < ad_model.pA.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.pA[i].data(), ed_model.pA[i].data(), 6, PREC );
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.pA.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.pA[i][idir].data(), ed_d_model.pA[i].col(idir).data(), 6, 5e-5);
    }
  }

  // nominal check
  CHECK_CLOSE(ad_model.u.rows(), ed_model.u.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_model.u.data(), ed_model.u.data(), ad_model.u.rows(), PREC);
  // derivative check
  CHECK_ARRAY_CLOSE(ad_d_model.u.data(), ed_d_model.u.data(), ad_d_model.u.cols() * ad_model.u.rows(), 5e-5);

  // nominal check
  CHECK_CLOSE(ad_model.f.size(), ed_model.f.size(), PREC);
  for (unsigned i = 0; i < ad_model.f.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_model.f[i].data(), ed_model.f[i].data(), 6, PREC);
  }
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    for (unsigned i = 0; i < ad_model.f.size(); i++) {
      CHECK_ARRAY_CLOSE(ad_d_model.f[i][idir].data(), ed_d_model.f[i].col(idir).data(), 6, 5e-5);
    }
  }
*/
  return;
}


void checkConstraintSetsADvsFD (
    unsigned ndirs,
    ConstraintSet const & ad_cs,
    ADConstraintSet const & ad_d_cs,
    ConstraintSet const & fd_cs,
    ADConstraintSet const & fd_d_cs,
    double const& PREC
) {

  // nominal check
  CHECK_CLOSE(ad_cs.H.rows(), fd_cs.H.rows(), PREC);
  CHECK_CLOSE(ad_cs.H.cols(), fd_cs.H.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.H.data(), fd_cs.H.data(),
                    fd_cs.H.rows() * fd_cs.H.cols(),
                    PREC);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.H[idir].data(), fd_d_cs.H[idir].data(),
                      fd_cs.H.rows() * fd_cs.H.cols(),
                      PREC);
  }

  // nominal check
  CHECK_CLOSE(ad_cs.G.rows(), fd_cs.G.rows(), PREC);
  CHECK_CLOSE(ad_cs.G.cols(), fd_cs.G.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.G.data(), fd_cs.G.data(),
                    fd_cs.G.rows() * fd_cs.G.cols(),
                    PREC);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.G[idir].data(), fd_d_cs.G[idir].data(),
                      fd_cs.G.rows() * fd_cs.G.cols(),
                      PREC);
  }

  // nominal check
  CHECK_CLOSE(ad_cs.A.rows(), fd_cs.A.rows(), PREC);
  CHECK_CLOSE(ad_cs.A.cols(), fd_cs.A.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.A.data(), fd_cs.A.data(),
                    fd_cs.A.rows() * fd_cs.A.cols(),
                    PREC);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.A[idir].data(), fd_d_cs.A[idir].data(),
                      fd_cs.A.rows() * fd_cs.A.cols(),
                      PREC);
  }

  // nominal check
  CHECK_CLOSE(ad_cs.Gi.rows(), fd_cs.Gi.rows(), PREC);
  CHECK_CLOSE(ad_cs.Gi.cols(), fd_cs.Gi.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.Gi.data(), fd_cs.Gi.data(),
                    fd_cs.Gi.rows() * fd_cs.Gi.cols(),
                    PREC);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.Gi[idir].data(), fd_d_cs.Gi[idir].data(),
                      fd_cs.Gi.rows() * fd_cs.Gi.cols(),
                      PREC);
  }

  // nominal check
  CHECK_CLOSE(ad_cs.GSpi.rows(), fd_cs.GSpi.rows(), PREC);
  CHECK_CLOSE(ad_cs.GSpi.cols(), fd_cs.GSpi.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.GSpi.data(), fd_cs.GSpi.data(),
                    fd_cs.GSpi.rows() * fd_cs.GSpi.cols(),
                    PREC);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.GSpi[idir].data(), fd_d_cs.GSpi[idir].data(),
                      fd_cs.GSpi.rows() * fd_cs.GSpi.cols(),
                      PREC);
  }

  // nominal check
  CHECK_CLOSE(ad_cs.GSsi.rows(), fd_cs.GSsi.rows(), PREC);
  CHECK_CLOSE(ad_cs.GSsi.cols(), fd_cs.GSsi.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.GSsi.data(), fd_cs.GSsi.data(),
                    fd_cs.GSsi.rows() * fd_cs.GSsi.cols(),
                    PREC);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.GSsi[idir].data(), fd_d_cs.GSsi[idir].data(),
                      fd_cs.GSsi.rows() * fd_cs.GSsi.cols(),
                      PREC);
  }

  // nominal check
  CHECK_CLOSE(ad_cs.C.rows(), fd_cs.C.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.C.data(), fd_cs.C.data(), ad_cs.C.rows(), PREC);
  // derivative check
  CHECK_CLOSE(ad_d_cs.C.rows(), fd_d_cs.C.rows(), PREC);
  CHECK_CLOSE(ad_d_cs.C.cols(), fd_d_cs.C.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_d_cs.C.data(), fd_d_cs.C.data(),
                    ad_d_cs.C.rows() * ad_d_cs.C.cols(), PREC);

  // nominal check
  CHECK_CLOSE(ad_cs.err.rows(), fd_cs.err.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.err.data(), fd_cs.err.data(), ad_cs.err.rows(), PREC);
  // derivative check
  CHECK_CLOSE(ad_d_cs.err.rows(), fd_d_cs.err.rows(), PREC);
  CHECK_CLOSE(ad_d_cs.err.cols(), fd_d_cs.err.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_d_cs.err.data(), fd_d_cs.err.data(),
                    ad_d_cs.err.rows() * ad_d_cs.err.cols(), PREC);

  // nominal check
  CHECK_CLOSE(ad_cs.errd.rows(), fd_cs.errd.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.errd.data(), fd_cs.errd.data(),
                    ad_cs.errd.rows(), PREC);
  // derivative check
  CHECK_CLOSE(ad_d_cs.errd.rows(), fd_d_cs.errd.rows(), PREC);
  CHECK_CLOSE(ad_d_cs.errd.cols(), fd_d_cs.errd.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_d_cs.errd.data(), fd_d_cs.errd.data(),
                    ad_d_cs.errd.rows() * ad_d_cs.errd.cols(), PREC);

  // nominal check
  CHECK_CLOSE(ad_cs.QDDot_0.rows(), fd_cs.QDDot_0.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.QDDot_0.data(), fd_cs.QDDot_0.data(),
                    ad_cs.QDDot_0.rows(), PREC);
  // derivative check
  CHECK_CLOSE(ad_d_cs.QDDot_0.rows(), fd_d_cs.QDDot_0.rows(), PREC);
  CHECK_CLOSE(ad_d_cs.QDDot_0.cols(), fd_d_cs.QDDot_0.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_d_cs.QDDot_0.data(), fd_d_cs.QDDot_0.data(),
                    ad_d_cs.QDDot_0.rows() * ad_d_cs.QDDot_0.cols(), PREC);

  // nominal check
  CHECK_CLOSE(ad_cs.gamma.rows(), fd_cs.gamma.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.gamma.data(), fd_cs.gamma.data(),
                    ad_cs.gamma.rows(), PREC);
  // derivative check
  CHECK_CLOSE(ad_d_cs.gamma.rows(), fd_d_cs.gamma.rows(), PREC);
  CHECK_CLOSE(ad_d_cs.gamma.cols(), fd_d_cs.gamma.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_d_cs.gamma.data(), fd_d_cs.gamma.data(),
                    ad_d_cs.gamma.rows() * ad_d_cs.gamma.cols(), PREC);

  // nominal check
  CHECK_CLOSE(ad_cs.d_u.rows(), fd_cs.d_u.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.d_u.data(), fd_cs.d_u.data(),
                    ad_cs.d_u.rows(), PREC);
  // derivative check
  CHECK_CLOSE(ad_d_cs.d_u.rows(), fd_d_cs.d_u.rows(), PREC);
  CHECK_CLOSE(ad_d_cs.d_u.cols(), fd_d_cs.d_u.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_d_cs.d_u.data(), fd_d_cs.d_u.data(),
                    ad_d_cs.d_u.rows() * ad_d_cs.d_u.cols(), PREC);

  // nominal check
  CHECK_CLOSE(ad_cs.d_a.size(), fd_cs.d_a.size(), PREC);
  for (unsigned i = 0; i < ad_cs.d_a.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_cs.d_a[i].data(), fd_cs.d_a[i].data(),
                      ad_cs.d_a[i].rows(), PREC);
    // derivative check
    CHECK_CLOSE(ad_d_cs.d_a[i].rows(), fd_d_cs.d_a[i].rows(), PREC);
    CHECK_CLOSE(ad_d_cs.d_a[i].cols(), fd_d_cs.d_a[i].cols(), PREC);
    CHECK_ARRAY_CLOSE(ad_d_cs.d_a[i].data(), fd_d_cs.d_a[i].data(),
                      ad_d_cs.d_a[i].rows() * ad_d_cs.d_a[i].cols(), PREC);
  }


}


void checkConstraintSetsADvsED (
    unsigned ndirs,
    ConstraintSet const & ad_cs,
    ADConstraintSet const & ad_d_cs,
    ConstraintSet const & ed_cs,
    EDConstraintSet const & ed_d_cs,
    double const& PREC
) {

  // nominal check
  CHECK_CLOSE(ad_cs.H.rows(), ed_cs.H.rows(), PREC);
  CHECK_CLOSE(ad_cs.H.cols(), ed_cs.H.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.H.data(), ed_cs.H.data(),
                    ed_cs.H.rows() * ed_cs.H.cols(),
                    PREC);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.H[idir].data(), ed_d_cs.H[idir].data(),
                      ed_cs.H.rows() * ed_cs.H.cols(),
                      PREC);
  }

  // nominal check
  CHECK_CLOSE(ad_cs.G.rows(), ed_cs.G.rows(), PREC);
  CHECK_CLOSE(ad_cs.G.cols(), ed_cs.G.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.G.data(), ed_cs.G.data(),
                    ed_cs.G.rows() * ed_cs.G.cols(),
                    PREC);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.G[idir].data(), ed_d_cs.G[idir].data(),
                      ed_cs.G.rows() * ed_cs.G.cols(),
                      PREC);
  }

  /*
  // nominal check
  CHECK_CLOSE(ad_cs.A.rows(), ed_cs.A.rows(), PREC);
  CHECK_CLOSE(ad_cs.A.cols(), ed_cs.A.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.A.data(), ed_cs.A.data(),
                    ed_cs.A.rows() * ed_cs.A.cols(),
                    PREC);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.A[idir].data(), ed_d_cs.A[idir].data(),
                      ed_cs.A.rows() * ed_cs.A.cols(),
                      PREC);
  }

  // nominal check
  CHECK_CLOSE(ad_cs.Gi.rows(), ed_cs.Gi.rows(), PREC);
  CHECK_CLOSE(ad_cs.Gi.cols(), ed_cs.Gi.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.Gi.data(), ed_cs.Gi.data(),
                    ed_cs.Gi.rows() * ed_cs.Gi.cols(),
                    PREC);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.Gi[idir].data(), ed_d_cs.Gi[idir].data(),
                      ed_cs.Gi.rows() * ed_cs.Gi.cols(),
                      PREC);
  }

  // nominal check
  CHECK_CLOSE(ad_cs.GSpi.rows(), ed_cs.GSpi.rows(), PREC);
  CHECK_CLOSE(ad_cs.GSpi.cols(), ed_cs.GSpi.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.GSpi.data(), ed_cs.GSpi.data(),
                    ed_cs.GSpi.rows() * ed_cs.GSpi.cols(),
                    PREC);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.GSpi[idir].data(), ed_d_cs.GSpi[idir].data(),
                      ed_cs.GSpi.rows() * ed_cs.GSpi.cols(),
                      PREC);
  }

  // nominal check
  CHECK_CLOSE(ad_cs.GSsi.rows(), ed_cs.GSsi.rows(), PREC);
  CHECK_CLOSE(ad_cs.GSsi.cols(), ed_cs.GSsi.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.GSsi.data(), ed_cs.GSsi.data(),
                    ed_cs.GSsi.rows() * ed_cs.GSsi.cols(),
                    PREC);
  // derivative check
  for (unsigned idir = 0; idir < ndirs; idir++) {
    CHECK_ARRAY_CLOSE(ad_d_cs.GSsi[idir].data(), ed_d_cs.GSsi[idir].data(),
                      ed_cs.GSsi.rows() * ed_cs.GSsi.cols(),
                      PREC);
  }
*/

  // nominal check
  CHECK_CLOSE(ad_cs.C.rows(), ed_cs.C.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.C.data(), ed_cs.C.data(), ad_cs.C.rows(), PREC);
  // derivative check
  CHECK_CLOSE(ad_d_cs.C.rows(), ed_d_cs.C.rows(), PREC);
  CHECK_CLOSE(ad_d_cs.C.cols(), ed_d_cs.C.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_d_cs.C.data(), ed_d_cs.C.data(),
                    ad_d_cs.C.rows() * ad_d_cs.C.cols(), PREC);

/*
  // nominal check
  CHECK_CLOSE(ad_cs.err.rows(), ed_cs.err.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.err.data(), ed_cs.err.data(), ad_cs.err.rows(), PREC);
  // derivative check
  CHECK_CLOSE(ad_d_cs.err.rows(), ed_d_cs.err.rows(), PREC);
  CHECK_CLOSE(ad_d_cs.err.cols(), ed_d_cs.err.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_d_cs.err.data(), ed_d_cs.err.data(),
                    ad_d_cs.err.rows() * ad_d_cs.err.cols(), PREC);

  // nominal check
  CHECK_CLOSE(ad_cs.errd.rows(), ed_cs.errd.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.errd.data(), ed_cs.errd.data(),
                    ad_cs.errd.rows(), PREC);
  // derivative check
  CHECK_CLOSE(ad_d_cs.errd.rows(), ed_d_cs.errd.rows(), PREC);
  CHECK_CLOSE(ad_d_cs.errd.cols(), ed_d_cs.errd.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_d_cs.errd.data(), ed_d_cs.errd.data(),
                    ad_d_cs.errd.rows() * ad_d_cs.errd.cols(), PREC);

  // nominal check
  CHECK_CLOSE(ad_cs.QDDot_0.rows(), ed_cs.QDDot_0.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.QDDot_0.data(), ed_cs.QDDot_0.data(),
                    ad_cs.QDDot_0.rows(), PREC);
  // derivative check
  CHECK_CLOSE(ad_d_cs.QDDot_0.rows(), ed_d_cs.QDDot_0.rows(), PREC);
  CHECK_CLOSE(ad_d_cs.QDDot_0.cols(), ed_d_cs.QDDot_0.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_d_cs.QDDot_0.data(), ed_d_cs.QDDot_0.data(),
                    ad_d_cs.QDDot_0.rows() * ad_d_cs.QDDot_0.cols(), PREC);

  */
  // nominal check
  CHECK_CLOSE(ad_cs.gamma.rows(), ed_cs.gamma.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.gamma.data(), ed_cs.gamma.data(),
                    ad_cs.gamma.rows(), PREC);
  // derivative check
  CHECK_CLOSE(ad_d_cs.gamma.rows(), ed_d_cs.gamma.rows(), PREC);
  CHECK_CLOSE(ad_d_cs.gamma.cols(), ed_d_cs.gamma.cols(), PREC);
  RigidBodyDynamics::Math::MatrixNd v_error
    = (ad_d_cs.gamma - ed_d_cs.gamma).cwiseAbs();
  const double error = v_error.maxCoeff();
  if (error > PREC) {
    std::cout << "error         = " << error << std::endl;
    std::cout << "gamma_error   = " << std::endl << v_error << std::endl;
    std::cout << "ad_d_cs.gamma = " << std::endl << ad_d_cs.gamma << std::endl;
    std::cout << "ed_d_cs.gamma = " << std::endl << ed_d_cs.gamma << std::endl;
    std::cout << endl;
  }
  CHECK_ARRAY_CLOSE(ad_d_cs.gamma.data(), ed_d_cs.gamma.data(),
                    ad_d_cs.gamma.rows() * ad_d_cs.gamma.cols(), PREC);

  /*
  // nominal check
  CHECK_CLOSE(ad_cs.d_u.rows(), ed_cs.d_u.rows(), PREC);
  CHECK_ARRAY_CLOSE(ad_cs.d_u.data(), ed_cs.d_u.data(),
                    ad_cs.d_u.rows(), PREC);
  // derivative check
  CHECK_CLOSE(ad_d_cs.d_u.rows(), ed_d_cs.d_u.rows(), PREC);
  CHECK_CLOSE(ad_d_cs.d_u.cols(), ed_d_cs.d_u.cols(), PREC);
  CHECK_ARRAY_CLOSE(ad_d_cs.d_u.data(), ed_d_cs.d_u.data(),
                    ad_d_cs.d_u.rows() * ad_d_cs.d_u.cols(), PREC);

  // nominal check
  CHECK_CLOSE(ad_cs.d_a.size(), ed_cs.d_a.size(), PREC);
  for (unsigned i = 0; i < ad_cs.d_a.size(); i++) {
    CHECK_ARRAY_CLOSE(ad_cs.d_a[i].data(), ed_cs.d_a[i].data(),
                      ad_cs.d_a[i].rows(), PREC);
    // derivative check
    CHECK_CLOSE(ad_d_cs.d_a[i].rows(), ed_d_cs.d_a[i].rows(), PREC);
    CHECK_CLOSE(ad_d_cs.d_a[i].cols(), ed_d_cs.d_a[i].cols(), PREC);
    CHECK_ARRAY_CLOSE(ad_d_cs.d_a[i].data(), ed_d_cs.d_a[i].data(),
                      ad_d_cs.d_a[i].rows() * ad_d_cs.d_a[i].cols(), PREC);
  }

*/

}

// -----------------------------------------------------------------------------
} //namespace RigidBodyDynamics
// -----------------------------------------------------------------------------


