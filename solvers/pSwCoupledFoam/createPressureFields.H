// Assume this is total flux
#include "createPhi.H"

// Pressure flux
surfaceScalarField phiwP = ("phiwP", -Mwf*fvc::snGrad(p) * mesh.magSf());
surfaceScalarField phinP = ("phinP", -Mnf*fvc::snGrad(p) * mesh.magSf());

// Gravitational flux
surfaceScalarField phiwG("phiG", (Lwf * g) & mesh.Sf());
surfaceScalarField phinG("phiG", (Lnf * g) & mesh.Sf());

// Equation coeffs
volScalarField wStorage("wStorage", porosity*rFVFw);
volScalarField oStorage("oStorage", porosity*rFVFn);
