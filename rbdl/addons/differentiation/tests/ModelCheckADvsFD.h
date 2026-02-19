#ifndef MODELCHECKADVSFD_h
#define MODELCHECKADVSFD_h

#include <rbdl/rbdl.h>
#include "ModelAD.h"
#include "ModelED.h"
#include "ConstraintsAD.h"
#include "ConstraintsED.h"

// -----------------------------------------------------------------------------
namespace RigidBodyDynamics {
// -----------------------------------------------------------------------------

void checkModelsADvsFD (
    unsigned ndirs,
    Model const & ad_model,
    ADModel const & ad_d_model,
    Model const & fd_model,
    ADModel const & fd_d_model,
    double const & PREC = 1e-6
);

void checkModelsADvsED (
    unsigned ndirs,
    Model const & ad_model,
    ADModel const & ad_d_model,
    Model const & ed_model,
    EDModel const & ed_d_model,
    double const & PREC = 1e-12);

void checkConstraintSetsADvsFD (
    unsigned ndirs,
    ConstraintSet const & ad_cs,
    ADConstraintSet const & ad_d_cs,
    ConstraintSet const & fd_cs,
    ADConstraintSet const & fd_d_cs,
    double const & PREC = 1e-6);

void checkConstraintSetsADvsED (
    unsigned ndirs,
    ConstraintSet const & ad_cs,
    ADConstraintSet const & ad_d_cs,
    ConstraintSet const & fd_cs,
    EDConstraintSet const & fd_d_cs,
    double const & PREC = 1e-12);

// -----------------------------------------------------------------------------
} // namespace RigidBodyDynamics
// -----------------------------------------------------------------------------

#endif // MODELCHECKADVSFD_h
