// Access to relative permeability
krModel->correct();
const volScalarField& krn = krModel->phase2Kr();
const volScalarField& krw = krModel->phase1Kr();
const volScalarField& dkrndS = krModel->phase2dKrdS();
const volScalarField& dkrwdS = krModel->phase1dKrdS();

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
surfaceTensorField Mnf ("Mnf",Kf*krnf/munf);
surfaceTensorField Lnf ("Lnf",rhonf*Kf*krnf/munf);	
surfaceTensorField Mwf ("Mwf",Kf*krwf/muwf);
surfaceTensorField Lwf ("Lwf",rhowf*Kf*krwf/muwf);
surfaceTensorField Mf ("Mf",Mnf+Mwf);
surfaceTensorField Lf ("Lf",Lnf+Lwf);
//surfaceTensorField Fwf ("Fwf",Mwf/Mf); There is no / operator between tensors, hhh
surfaceScalarField Fwf ("Fwf",(krwf/muwf) / ((krnf/munf) + (krwf/muwf)));
volScalarField Fw ("Fw",(krw/muw) / ( (krn/mun) + (krw/muw) ));

// Capillary flux
pcModel->correct();
// snGrad  is out of place
//surfaceScalarField phiPc("phiPc",Mwf*fvc::interpolate(pcModel->dpcdS(),"pc")*fvc::snGrad(phasew.alpha())*mesh.magSf());
// so, the grad is calculated at cells then interpolated 
// I hope this doesn't generate bugs, take a deeper look at snGrad
surfaceVectorField gradPc("gradPc",fvc::interpolate(pcModel->dpcdS()*fvc::grad(phasew.alpha()),"pc"));
surfaceScalarField phiPc("phiPc",(Mwf & gradPc) & mesh.Sf());
