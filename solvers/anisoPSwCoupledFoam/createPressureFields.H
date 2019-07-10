// Assume this is total flux
#include "createPhi.H"

// Pressure flux
// Need to interpolate volume grad to faces, is there a better way??
surfaceVectorField gradP = ("gradP", fvc::interpolate(fvc::grad(p),"gradP"));
surfaceScalarField phiwP = ("phiwP", (-Mwf & gradP) & mesh.Sf());
surfaceScalarField phinP = ("phinP", (-Mnf & gradP) & mesh.Sf());
/////surfaceScalarField phiwP = ("phiwP", -Mwf * fvc::snGrad(p) * mesh.magSf());
/////surfaceScalarField phinP = ("phinP", -Mnf * fvc::snGrad(p) * mesh.magSf());

// Gravitational flux
surfaceScalarField phiwG("phiG", (Lwf & g) & mesh.Sf());
surfaceScalarField phinG("phiG", (Lnf & g) & mesh.Sf());

// Equation coeffs
volScalarField wStorage("wStorage", porosity*rFVFw);
volScalarField oStorage("oStorage", porosity*rFVFn);
