// Access to relative permeability
krModel->correct();
volScalarField krc
(
    canPhase.name()+".kr",
    krModel->kr(krModel->phaseIndex(canPhase.name()))
);
volScalarField krn
(
    noncanPhase.name()+".kr",
    krModel->kr(krModel->phaseIndex(noncanPhase.name()))
);
volScalarField dkrcdS
(
    canPhase.name()+".dkrdS",
    krModel->dkrdS(krModel->phaseIndex(canPhase.name()))
);
volScalarField dkrndS
(
    noncanPhase.name()+".dkrdS",
    krModel->dkrdS(krModel->phaseIndex(noncanPhase.name()))
);

// Interpolations to faces, profile this
surfaceScalarField krcf ("krcf", fvc::interpolate(krc, krc.name()));
surfaceScalarField krnf ("krnf", fvc::interpolate(krn, krn.name()));

surfaceScalarField dkrcfdS ("dkrcfdS", fvc::interpolate(dkrcdS, dkrcdS.name()));
surfaceScalarField dkrnfdS ("dkrnfdS", fvc::interpolate(dkrndS, dkrndS.name()));

// Mobilities and flow fractions
surfaceScalarField Mcf ("Mcf", Kf*krcf/muc);
surfaceScalarField Lcf ("Lcf", rhoc*Kf*krcf/muc);
surfaceScalarField Mnf ("Mnf", Kf*krnf/mun);
surfaceScalarField Lnf ("Lnf", rhon*Kf*krnf/mun);
surfaceScalarField Mf ("Mf", Mnf+Mcf);
surfaceScalarField Lf ("Lf", Lnf+Lcf);
surfaceScalarField Fcf ("Fcf", Mcf/Mf);
volScalarField Fc ("Fc", (krc/muc) / ( (krn/mun) + (krc/muc) ));

// Capillary flux
pcModel->correct();
surfaceScalarField phiPc
(
    "phiPc",
    Mcf * fvc::interpolate(pcModel->dpcdS(),"pc")
    * fvc::snGrad(canPhase.alpha())*mesh.magSf()
);
