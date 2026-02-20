#ifndef FDMODELENTRY_FDC_H
#define FDMODELENTRY_FDC_H

#include "rbdl/rbdl.h"

#include "ModelAD.h"
#include "ConstraintsAD.h"

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
);

void computeFDCEntry(
  ConstraintSet const &csph,
  ConstraintSet const &csmh,
  const double &H,
  const int &idir,
  ADConstraintSet &fd_cs
);

// -----------------------------------------------------------------------------
} // FDC
// -----------------------------------------------------------------------------
} // RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif // FDMODELENTRY_FDC_H
