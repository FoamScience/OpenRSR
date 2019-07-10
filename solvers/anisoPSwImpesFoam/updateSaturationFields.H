// Update relative permeabilities
krModel->correct();

// Interpolate kr data
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
//Fwf = Mwf/(Mf+VSMALLF);
Fwf = (krwf/muwf) / ((krnf/munf) + (krwf/muwf));

// Fractional flow coeff.
Fw = (krw/muw) / ( (krn/mun) + (krw/muw) );

// Capillary Flux
pcModel->correct();
gradPc = fvc::interpolate(pcModel->dpcdS()*fvc::grad(phasew.alpha()),"pc");
phiPc = (Mwf & gradPc) & mesh.Sf();
