// Update fluxes
phiG = (Lf * g) & mesh.Sf();
phiwG = (Lwf*g) & mesh.Sf();

// Reconstruct total flux
// phiP is calculated in pEqn.H
phi = phiP+phiG+phiPc;

// Continuity errors
#include "continuityErrs.H"

// Overwrite BCs
phiw == Fwf*phiP + phiwG + phiPc;
phin == phi - phiw;

// Reconstruct velocity fields
U = fvc::reconstruct(phi);
U.correctBoundaryConditions();

Uw = fvc::reconstruct(phiw);
Un = U-Uw;

Uw.correctBoundaryConditions();  
Un.correctBoundaryConditions();


