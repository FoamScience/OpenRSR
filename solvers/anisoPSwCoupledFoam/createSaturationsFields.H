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
///// Transmissibility scalarFields to tensorFields
surfaceTensorField Mnf ("Mnf",Kf*krnf*rFVFnf/munf);
surfaceTensorField Lnf ("Lnf",rhonf*Mnf);	
surfaceTensorField Mwf ("Mwf",Kf*krwf*rFVFwf/muwf);
surfaceTensorField Lwf ("Lwf",rhowf*Mwf);
surfaceTensorField Mf ("Mf",Mnf+Mwf);
surfaceTensorField Lf ("Lf",Lnf+Lwf);
////surfaceTensorField Fwf ("Fwf",Mwf/(Mf+VSMALLF));
surfaceScalarField Fwf ("Fwf", (krwf*rFVFwf/muwf) / ((krnf*rFVFnf/munf) + (krwf*rFVFwf/muwf)));
volScalarField Fw ("Fw",(krw*rFVFw/muw) / ( (krn*rFVFn/mun) + (krw*rFVFw/muw) ));

// Capillary flux
pcModel->correct();
surfaceVectorField gradPc("gradPc",fvc::interpolate(pcModel->dpcdS()*fvc::grad(phasew.alpha()),"pc"));
surfaceScalarField phiPc("phiPc",(Mwf & gradPc) & mesh.Sf());
/////surfaceScalarField phiPc("phiPc",Mwf*fvc::interpolate(pcModel->dpcdS(),"pc")*fvc::snGrad(phasew.alpha())*mesh.magSf());
