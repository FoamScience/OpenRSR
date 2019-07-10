// Update relative permeabilities
krModel->correct();
// Update phase stuff
mixture->correct();
// Update porosity
porosity = porosity_0*Foam::pow( mathematicalConstant::e, compressibility*(p-p0));

// Interpolate kr data
//krn = krModel->phase2Kr();
//krw = krModel->phase1Kr();
//dkrndS = krModel->phase2dKrdS();
//dkrwdS = krModel->phase1dKrdS();
krnf = fvc::interpolate(krn,"krn");
krwf = fvc::interpolate(krw,"krw");
dkrnfdS = fvc::interpolate(dkrndS,"krn");
dkrwfdS = fvc::interpolate(dkrwdS,"krw");
rhonf = fvc::interpolate(rhon,"rho");
rhowf = fvc::interpolate(rhow,"rho");
rFVFnf = fvc::interpolate(rFVFn,"rFVF");
rFVFwf = fvc::interpolate(rFVFw,"rFVF");

// Update mobility and Transmissiblity Data
Mnf = Kf*krnf*rFVFnf/munf;
Lnf = rhonf*Mnf;
Mwf = Kf*krwf*rFVFwf/muwf;
Lwf = rhowf*Mwf;
Mf = Mnf+Mwf;
Lf = Lnf+Lwf;
/////Fwf = Mwf/(Mf+VSMALLF);
Fwf = (krwf*rFVFwf/muwf) / ((krnf*rFVFnf/munf) + (krwf*rFVFwf/muwf));
Fw = (krw*rFVFw/muw) / ( (krn*rFVFn/mun) + (krw*rFVFw/muw) );

// Capillary Flux
pcModel->correct();
/////phiPc = Mwf * fvc::interpolate(pcModel->dpcdS(),"pc")* fvc::snGrad(phasew.alpha()) * mesh.magSf();
gradPc = fvc::interpolate(pcModel->dpcdS()*fvc::grad(phasew.alpha()),"pc");
phiPc = (Mwf & gradPc) & mesh.Sf();
