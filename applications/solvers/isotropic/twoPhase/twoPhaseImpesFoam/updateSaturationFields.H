// Update non-canonical saturation
noncanPhase.alpha() = 1 - canPhase.alpha();

// Update relative permeabilities
krModel->correct();

// Re-Interpolate kr data
krc = krModel->kr(krModel->phaseIndex(canPhase.name()));
krn = krModel->kr(krModel->phaseIndex(noncanPhase.name()));
dkrcdS = krModel->dkrdS(krModel->phaseIndex(canPhase.name()));
dkrndS = krModel->dkrdS(krModel->phaseIndex(noncanPhase.name()));

krcf = fvc::interpolate(krc, canPhase.name()+".kr");
krnf = fvc::interpolate(krn, noncanPhase.name()+".kr");
dkrcfdS = fvc::interpolate(dkrcdS, canPhase.name()+".dkrdS");
dkrnfdS = fvc::interpolate(dkrndS, noncanPhase.name()+".dkrdS");

// Update mobility and Transmissiblity Data
Mnf = Kf*krnf/mun;
Lnf = rhon*Mnf;	
Mcf = Kf*krcf/muc;
Lcf = rhoc*Mcf;
Mf = Mnf+Mcf;
Lf = Lnf+Lcf;
Fcf = Mcf/Mf;

// Fractional flow coeff.
Fc = (krc/muc) / ( (krc/muc) + (krc/muc) );

// Capillary Flux
pcModel->correct();
phiPc = Mcf * fvc::interpolate(pcModel->dpcdS(),"pc")
        * fvc::snGrad(canPhase.alpha()) * mesh.magSf();
