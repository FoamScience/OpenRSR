// Update gravitation fluxes
phiG = (Lf * g) & mesh.Sf();
phicG = (Lcf*g) & mesh.Sf();

// Reconstruct total flux
// phiP is calculated in pEqn.H
phi == phiP+phiG+phiPc;

// Continuity errors (on global flux: phi)
#include "continuityErrs.H"

// Overwrite BCs
phic == Fcf*phiP + phicG + phiPc;
phin == phi - phic;

// Reconstruct velocity fields
U = fvc::reconstruct(phi);
U.correctBoundaryConditions();

Uc = fvc::reconstruct(phic);
Un = U-Uc;

Uc.correctBoundaryConditions();
Un.correctBoundaryConditions();
