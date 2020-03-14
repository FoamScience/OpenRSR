// FROM porousMultiphaseFoam project by P. Horgue et al.
// Adjust to suit the solver and fixed some bugs in aniso version
{
    // Inertia
    volScalarField dfc (
        "dfc",
        (dkrcdS*krn - dkrndS*krc)
        // Should this be stored first?
        /(muc/mun*Foam::pow(krn,2) + 2*krn*krc + mun/muc*Foam::pow(krc,2))
    );

    dimensionedScalar smallRate("smallRate",dimVolume/dimTime, SMALL);

    dfc -= 
        K*(rhon-rhoc)
        *fvc::surfaceSum(mag(mesh.Sf() & g))
        /fvc::surfaceSum(mag(phi)+smallRate)
        *(Foam::pow(krn,2)*dkrcdS/mun + Foam::pow(krc,2)*dkrndS/muc)
        /(muc/mun*Foam::pow(krn,2) + 2*krn*krc + mun/muc*Foam::pow(krc,2));

    scalarField CFLCoats
    (
        runTime.deltaT()*dfc*fvc::surfaceSum(mag(phi))/porosity/mesh.V()
    );

    // Capillarity
    CFLCoats += 
        (runTime.deltaT()/porosity/mesh.V())
        *2*mag(pcModel->dpcdS())
        *fvc::surfaceSum(Kf*mesh.magSf()/mag(delta))
        *(krn*krc/(muc*krn+mun*krc));

    Info<< "Coats Number mean: " << gAverage(CFLCoats) << " max: " << gMax(CFLCoats) << endl;

    CFLUse = gMax(CFLCoats);
    maxDeltaTFact = maxCo/(CFLUse + SMALL);
      
}
