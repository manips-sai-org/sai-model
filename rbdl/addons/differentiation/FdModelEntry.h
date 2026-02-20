#ifndef FDMODELENTRY_h
#define FDMODELENTRY_h

#include "rbdl/rbdl.h"

#include "ModelAD.h"
#include "ConstraintsAD.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

void computeFDEntry(
    Model const & model,
    Model const & modelh,
    double h,
    int idir,
    ADModel & fd_model);

void computeFDEntry(
    ConstraintSet const & cs,
    ConstraintSet const & csh,
    double h, int idir,
    ADConstraintSet &fd_cs);

// -----------------------------------------------------------------------------
} // RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif // FDMODELENTRY_h
