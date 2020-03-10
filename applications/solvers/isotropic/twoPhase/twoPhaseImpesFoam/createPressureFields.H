// Assume this is total flux
Info<< nl;
#include "createPhi.H"

// Flux generated from pressure difference
surfaceScalarField phiP = phi;

// Gravitational fluxes
surfaceScalarField phiG("phiG", (Lf * g) & mesh.Sf());
surfaceScalarField phicG(canPhase.name()+".phiG", (Lcf * g) & mesh.Sf());
