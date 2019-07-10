// For convenience
dimensionedScalar VSMALLF ("VSMALLF", dimensionSet(-1,3,1,0,0,0,0), VSMALL);

// Access to relative permeability
krModel->correct();
volScalarField krn = krModel->phase2Kr();
volScalarField krw = krModel->phase1Kr();
volScalarField dkrndS = krModel->phase2dKrdS();
volScalarField dkrwdS = krModel->phase1dKrdS();

// Interpolations to faces, profile this
surfaceScalarField krnf ("krnf",fvc::interpolate(krn,"krn"));
surfaceScalarField krwf ("krwf",fvc::interpolate(krw,"krw"));
surfaceScalarField dkrnfdS ("dkrnfdS",fvc::interpolate(dkrndS,"krn"));
surfaceScalarField dkrwfdS ("dkrwfdS",fvc::interpolate(dkrwdS,"krw"));
surfaceScalarField munf ("munf",fvc::interpolate(mun,"mu"));
surfaceScalarField muwf ("muwf",fvc::interpolate(muw,"mu"));
surfaceScalarField rhonf ("rhonf",fvc::interpolate(rhon,"rho"));
surfaceScalarField rhowf ("rhowf",fvc::interpolate(rhow,"rho"));

// Mobilities and flow fractions
surfaceScalarField Mnf ("Mnf",Kf*krnf/munf);
surfaceScalarField Lnf ("Lnf",rhonf*Kf*krnf/munf);	
surfaceScalarField Mwf ("Mwf",Kf*krwf/muwf);
surfaceScalarField Lwf ("Lwf",rhowf*Kf*krwf/muwf);
surfaceScalarField Mf ("Mf",Mnf+Mwf);
surfaceScalarField Lf ("Lf",Lnf+Lwf);
surfaceScalarField Fwf ("Fbf",Mwf/Mf);
volScalarField Fw ("Fw",(krw/muw) / ( (krn/mun) + (krw/muw) ));

// Capillary flux
pcModel->correct();
surfaceScalarField phiPc("phiPc",Mwf*fvc::interpolate(pcModel->dpcdS(),"pc")*fvc::snGrad(phasew.alpha())*mesh.magSf());
