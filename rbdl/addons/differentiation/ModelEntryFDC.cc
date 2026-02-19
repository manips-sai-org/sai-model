#include "ModelEntryFDC.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------
namespace FDC {
// -----------------------------------------------------------------------------

void computeFDCEntry(
  Model const &modelph,
  Model const &modelmh,
  const double &H,
  const int &idir,
  ADModel &fd_model
) {
  // directly compute right denominator
  const double H2 = 2*H;

  for (unsigned i = 0; i < modelmh.X_lambda.size(); i++) {
    fd_model.X_lambda[i][idir].E
      = (modelph.X_lambda[i].E - modelmh.X_lambda[i].E) / H2;

    fd_model.X_lambda[i][idir].r
      = (modelph.X_lambda[i].r - modelmh.X_lambda[i].r) / H2;
  }
  for (unsigned i = 0; i < modelmh.X_base.size(); i++) {
    fd_model.X_base[i][idir].E
      = (modelph.X_base[i].E - modelmh.X_base[i].E) / H2;

    fd_model.X_base[i][idir].r
      = (modelph.X_base[i].r - modelmh.X_base[i].r) / H2;
  }
  for (unsigned i = 0; i < modelmh.X_J.size(); i++) {
    fd_model.X_J[i][idir].E = (modelph.X_J[i].E - modelmh.X_J[i].E) / H2;
    fd_model.X_J[i][idir].r = (modelph.X_J[i].r - modelmh.X_J[i].r) / H2;
  }
  for (unsigned i = 0; i < modelmh.Ic.size(); i++) {
    fd_model.Ic[i][idir].m = (modelph.Ic[i].m - modelmh.Ic[i].m) / H2;
    fd_model.Ic[i][idir].h = (modelph.Ic[i].h - modelmh.Ic[i].h) / H2;
    fd_model.Ic[i][idir].Ixx = (modelph.Ic[i].Ixx - modelmh.Ic[i].Ixx) / H2;
    fd_model.Ic[i][idir].Iyx = (modelph.Ic[i].Iyx - modelmh.Ic[i].Iyx) / H2;
    fd_model.Ic[i][idir].Iyy = (modelph.Ic[i].Iyy - modelmh.Ic[i].Iyy) / H2;
    fd_model.Ic[i][idir].Izx = (modelph.Ic[i].Izx - modelmh.Ic[i].Izx) / H2;
    fd_model.Ic[i][idir].Izy = (modelph.Ic[i].Izy - modelmh.Ic[i].Izy) / H2;
    fd_model.Ic[i][idir].Izz = (modelph.Ic[i].Izz - modelmh.Ic[i].Izz) / H2;
  }
  for (unsigned i = 0; i < modelmh.IA.size(); i++) {
    fd_model.IA[i][idir] = (modelph.IA[i] - modelmh.IA[i]) / H2;
  }
  for (unsigned i = 0; i < modelmh.c.size(); i++) {
    fd_model.c[i][idir] = (modelph.c[i] - modelmh.c[i]) / H2;
  }
  for (unsigned i = 0; i < modelmh.c_J.size(); i++) {
    fd_model.c_J[i][idir] = (modelph.c_J[i] - modelmh.c_J[i]) / H2;
  }
  for (unsigned i = 0; i < modelmh.v.size(); i++) {
    fd_model.v[i][idir] = (modelph.v[i] - modelmh.v[i]) / H2;
  }
  for (unsigned i = 0; i < modelmh.v_J.size(); i++) {
    fd_model.v_J[i][idir] = (modelph.v_J[i] - modelmh.v_J[i]) / H2;
  }
  for (unsigned i = 0; i < modelmh.a.size(); i++) {
    fd_model.a[i][idir] = (modelph.a[i] - modelmh.a[i]) / H2;
  }
  for (unsigned i = 0; i < modelmh.U.size(); i++) {
    fd_model.U[i][idir] = (modelph.U[i] - modelmh.U[i]) / H2;
  }
  for (unsigned i = 0; i < modelmh.pA.size(); i++) {
    fd_model.pA[i][idir] = (modelph.pA[i] - modelmh.pA[i]) / H2;
  }
  for (unsigned i = 0; i < modelmh.u.rows(); i++) {
    fd_model.u(i, idir) = (modelph.u(i) - modelmh.u(i)) / H2;
  }
  for (unsigned i = 0; i < modelmh.d.rows(); i++) {
    fd_model.d(i, idir) = (modelph.d(i) - modelmh.d(i)) / H2;
  }
  for (unsigned i = 0; i < modelmh.f.size(); i++) {
    fd_model.f[i][idir] = (modelph.f[i] - modelmh.f[i]) / H2;
  }
}

void computeFDCEntry(
  ConstraintSet const &csph,
  ConstraintSet const &csmh,
  const double &H,
  const int &idir,
  ADConstraintSet &fd_cs
) {
  // directly compute right denominator
  const double H2 = 2*H;

  fd_cs.G[idir]  = (csph.G - csmh.G) / H2;
  fd_cs.Gi[idir] = (csph.Gi - csmh.Gi) / H2;
  fd_cs.GSpi[idir] = (csph.GSpi - csmh.GSpi) / H2;
  fd_cs.GSsi[idir] = (csph.GSsi - csmh.GSsi) / H2;
  fd_cs.A[idir]  = (csph.A - csmh.A) / H2;
  fd_cs.H[idir]  = (csph.H - csmh.H) / H2;
  fd_cs.b.col(idir) = (csph.b - csmh.b) / H2;
  fd_cs.v_plus.col(idir) = (csph.v_plus - csmh.v_plus) / H2;
  fd_cs.x.col(idir) = (csph.x - csmh.x) / H2;
  fd_cs.impulse.col(idir) = (csph.impulse - csmh.impulse) / H2;
  fd_cs.QDDot_0.col(idir) = (csph.QDDot_0 - csmh.QDDot_0) / H2;
  fd_cs.C.col(idir) = (csph.C - csmh.C) / H2;
  fd_cs.gamma.col(idir) = (csph.gamma - csmh.gamma) / H2;
  fd_cs.force.col(idir) = (csph.force - csmh.force) / H2;
  fd_cs.err.col(idir) = (csph.err - csmh.err) / H2;
  fd_cs.errd.col(idir) = (csph.errd - csmh.errd) / H2;
}

// -----------------------------------------------------------------------------
} // FDC
// -----------------------------------------------------------------------------
} // RigidBodyDynamics
// -----------------------------------------------------------------------------
