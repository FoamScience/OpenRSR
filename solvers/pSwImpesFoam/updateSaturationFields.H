// Update relative permeabilities
krModel->correct();

// Interpolate kr data
krn = krModel->phase2Kr();
krw = krModel->phase1Kr();
dkrndS = krModel->phase2dKrdS();
dkrwdS = krModel->phase1dKrdS();
krnf = fvc::interpolate(krn,"krn");
krwf = fvc::interpolate(krw,"krw");
dkrnfdS = fvc::interpolate(dkrndS,"krn");
dkrwfdS = fvc::interpolate(dkrwdS,"krw");

// Update mobility and Transmissiblity Data
Mnf = Kf*krnf/munf;
Lnf = rhonf*Kf*krnf/munf;	
Mwf = Kf*krwf/muwf;
Lwf = rhowf*Kf*krwf/muwf;
Mf = Mnf+Mwf;
Lf = Lnf+Lwf;
Fwf = Mwf/(Mf+VSMALLF);

// Fractional flow coeff.
Fw = (krw/muw) / ( (krn/mun) + (krw/muw) );

// Capillary Flux
pcModel->correct();
phiPc = Mwf * fvc::interpolate(pcModel->dpcdS(),"pc")* fvc::snGrad(phasew.alpha()) * mesh.magSf();

