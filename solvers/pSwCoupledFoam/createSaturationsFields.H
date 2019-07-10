// For convenience
dimensionedScalar VSMALLF ("VSMALLF", dimensionSet(-1,3,1,0,0,0,0), VSMALL);

// Access to relative permeability
krModel->correct();
volScalarField& krn = krModel->phase2Kr();
volScalarField& krw = krModel->phase1Kr();
volScalarField& dkrndS = krModel->phase2dKrdS();
volScalarField& dkrwdS = krModel->phase1dKrdS();

// Interpolations to faces, profile this
surfaceScalarField krnf ("krnf",fvc::interpolate(krn,"krn"));
surfaceScalarField krwf ("krwf",fvc::interpolate(krw,"krw"));
surfaceScalarField dkrnfdS ("dkrnfdS",fvc::interpolate(dkrndS,"krn"));
surfaceScalarField dkrwfdS ("dkrwfdS",fvc::interpolate(dkrwdS,"krw"));
surfaceScalarField munf ("munf",fvc::interpolate(mun,"mu"));
surfaceScalarField muwf ("muwf",fvc::interpolate(muw,"mu"));
surfaceScalarField rhonf ("rhonf",fvc::interpolate(rhon,"rho"));
surfaceScalarField rhowf ("rhowf",fvc::interpolate(rhow,"rho"));
surfaceScalarField rFVFnf ("rFVFnf",fvc::interpolate(rFVFn,"rFVF"));
surfaceScalarField rFVFwf ("rFVFwf",fvc::interpolate(rFVFw,"rFVF"));

// Mobilities and flow fractions
surfaceScalarField Mnf ("Mnf",Kf*krnf*rFVFnf/munf);
surfaceScalarField Lnf ("Lnf",rhonf*Mnf);	
surfaceScalarField Mwf ("Mwf",Kf*krwf*rFVFwf/muwf);
surfaceScalarField Lwf ("Lwf",rhowf*Mwf);
surfaceScalarField Mf ("Mf",Mnf+Mwf);
surfaceScalarField Lf ("Lf",Lnf+Lwf);
surfaceScalarField Fwf ("Fbf",Mwf/(Mf+VSMALLF));
volScalarField Fw ("Fw",(krw*rFVFw/muw) / ( (krn*rFVFn/mun) + (krw*rFVFw/muw) ));

// Capillary flux
pcModel->correct();
surfaceScalarField phiPc("phiPc",Mwf*fvc::interpolate(pcModel->dpcdS(),"pc")*fvc::snGrad(phasew.alpha())*mesh.magSf());

