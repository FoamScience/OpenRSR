// Assume this is total flux
#include "createPhi.H"

// Pressure flux
surfaceScalarField phiP = phi;

// Gravitational flux
surfaceScalarField phiG("phiG", (Lf * g) & mesh.Sf());
surfaceScalarField phicG("phiG", (Lcf * g) & mesh.Sf());
