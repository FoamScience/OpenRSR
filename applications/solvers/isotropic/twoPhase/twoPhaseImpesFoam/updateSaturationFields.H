noncanPhase.alpha() = 1 - canPhase.alpha();

// Update relative permeabilities
krModel->correct();

// Interpolate kr data
krc = krModel->kr(krModel->phaseIndex(canPhase.name()));
krn = krModel->kr(krModel->phaseIndex(noncanPhase.name()));
dkrcdS = krModel->dkrdS(krModel->phaseIndex(canPhase.name()));
dkrndS = krModel->dkrdS(krModel->phaseIndex(noncanPhase.name()));

krnf = fvc::interpolate(krn, krn.name());
krcf = fvc::interpolate(krc, krc.name());
dkrcfdS = fvc::interpolate(dkrcdS, dkrcdS.name());
dkrnfdS = fvc::interpolate(dkrndS, dkrndS.name());

// Update mobility and Transmissiblity Data
Mcf = Kf*krcf/muc;
Lcf = rhoc*Kf*krcf/muc;
Mnf = Kf*krnf/mun;
Lnf = rhon*Kf*krnf/mun;
Mf = Mnf+Mcf;
Lf = Lnf+Lcf;
Fcf = Mcf/Mf;

// Fractional flow coeff.
Fc = (krc/muc) / ( (krn/mun) + (krc/muc) );

// Capillary Flux
pcModel->correct();
phiPc = Mcf * fvc::interpolate(pcModel->dpcdS(),"pc")* fvc::snGrad(canPhase.alpha()) * mesh.magSf();

