// FROM porousMultiphaseFoam project by P. Horgue et al.
// Adjust to suit the solver and fixed some bugs in aniso version
{
    // Inertia
    volScalarField dfw (
        "dfw",
        (dkrwdS*krn - dkrndS*krw) 
        // Should this be stored first?
        /(muw/mun*Foam::pow(krn,2) + 2*krn*krw + mun/muw*Foam::pow(krw,2)) 
    );

    dimensionedScalar smallRate("smallRate",dimVolume/dimTime, SMALL);

    dfw -= 
        K*(rhon-rhow)
        *fvc::surfaceSum(mag(mesh.Sf() & g))
        /fvc::surfaceSum(mag(phi)+smallRate)
        *(Foam::pow(krn,2)*dkrwdS/mun + Foam::pow(krw,2)*dkrndS/muw)
        /(muw/mun*Foam::pow(krn,2) + 2*krn*krw + mun/muw*Foam::pow(krw,2));

    scalarField CFLCoats
    (
        runTime.deltaT()*dfw*fvc::surfaceSum(mag(phi))/porosity/mesh.V()
    );

    // Capillarity
    CFLCoats += 
        (runTime.deltaT()/porosity/mesh.V())
        *2*mag(pcModel->dpcdS())
        *fvc::surfaceSum(Kf*mesh.magSf()/mag(delta))
        *(krn*krw/(muw*krn+mun*krw));

    Info<< "Coats Number mean: " << gAverage(CFLCoats) << " max: " << gMax(CFLCoats) << endl;

    CFLUse = gMax(CFLCoats);
    maxDeltaTFact = maxCo/(CFLUse + SMALL);
      
}
