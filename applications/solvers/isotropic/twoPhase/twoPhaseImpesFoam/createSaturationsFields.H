// For convenience
dimensionedScalar VSMALLF ("VSMALLF", dimensionSet(-1,3,1,0,0,0,0), VSMALL);

// Check canonical phase consistency
if (canPhase.name() != krModel->canonicalPhase())
{
    WarningIn(__PRETTY_FUNCTION__)
        << "Relative Permeability Model has its canonical phase" 
        << " set to " << krModel->canonicalPhase()
        << " but the solver recognizes " << canPhase.name()
        << " as the canonical phase."
        << endl;
}

// Control Kr data
krModel->correct();
volScalarField krc = krModel->kr(krModel->phaseIndex(canPhase.name()));
volScalarField krn = krModel->kr(krModel->phaseIndex(noncanPhase.name()));
volScalarField dkrcdS = krModel->dkrdS(krModel->phaseIndex(canPhase.name()));
volScalarField dkrndS = krModel->dkrdS(krModel->phaseIndex(noncanPhase.name()));

// Interpolations to faces
surfaceScalarField krcf ("krcf", fvc::interpolate(krc, canPhase.name()+".kr"));
surfaceScalarField krnf ("krnf", fvc::interpolate(krn, noncanPhase.name()+".kr"));
surfaceScalarField dkrcfdS ("dkrcfdS", fvc::interpolate(dkrcdS, canPhase.name()+".dkrdS"));
surfaceScalarField dkrnfdS ("dkrnfdS", fvc::interpolate(dkrndS, noncanPhase.name()+".dkrdS"));

// Mobilities and flow fractions
surfaceScalarField Mnf ("Mnf",Kf*krnf/mun);
surfaceScalarField Lnf ("Lnf",rhon*Mnf);	
surfaceScalarField Mcf ("Mcf",Kf*krcf/muc);
surfaceScalarField Lcf ("Lcf",rhoc*Mcf);
surfaceScalarField Mf ("Mf",Mnf+Mcf);
surfaceScalarField Lf ("Lf",Lnf+Lcf);
surfaceScalarField Fcf ("Fcf",Mcf/Mf);

// Fractional flow
volScalarField Fc ("Fc",(krc/muc) / ( (krn/mun) + (krc/muc) ));

// Capillary flux
pcModel->correct();
surfaceScalarField phiPc
(
    "phiPc",
    Mcf * fvc::interpolate(pcModel->dpcdS(),"pc")
        * fvc::snGrad(canPhase.alpha())*mesh.magSf()
);
