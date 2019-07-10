// Update fluxes
//phiG = (Lf * g) & mesh.Sf();
phiwG = (Lwf & g) & mesh.Sf();
phinG = (Lnf & g) & mesh.Sf();

// Reconstruct total flux
// phiP is calculated in pEqn.H
/////phiwP = -Mwf*fvc::snGrad(p) * mesh.magSf();
/////phinP = -Mnf*fvc::snGrad(p) * mesh.magSf();
gradP = fvc::interpolate(fvc::grad(p),"gradP");
phiwP = (-Mwf & gradP) & mesh.Sf();
phinP = (-Mnf & gradP) & mesh.Sf();
phi == phiwP+phinP+phiwG+phinG+phiPc;

// Continuity errors
#include "continuityErrs.H"

// Overwrite BCs
phiw == phiwP + phiwG + phiPc;
phin == phinP + phinG;

// Reconstruct velocity fields
U = fvc::reconstruct(phi);
U.correctBoundaryConditions();

Uw = fvc::reconstruct(phiw);
Un = U-Uw;

Uw.correctBoundaryConditions();  
Un.correctBoundaryConditions();

// Equation coeffs
wStorage = porosity*rFVFw;
oStorage = porosity*rFVFn;
