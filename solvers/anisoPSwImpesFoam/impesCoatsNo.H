// FROM porousMultiphaseFoam project by P. Horgue et al.
{
    // Inertia
    volTensorField dfw (
        "dfw",
        tensor::I
        *(dkrwdS*krn - dkrndS*krw) 
        // Should this be stored first?
        /(muw/mun*Foam::pow(krn,2) + 2*krn*krw + mun/muw*Foam::pow(krw,2)) 
    );

    dimensionedScalar smallRate("smallRate",dimVolume/dimTime, SMALL);

    // dfw  is needed to remain scalar, so, use K at faces
    dfw -= 
        //K*(rhon-rhow)
        (rhon-rhow)
        //*fvc::surfaceSum(mag(mesh.Sf() & g))
        *fvc::surfaceSum(Kf * mag(mesh.Sf() & g))
        /fvc::surfaceSum(mag(phi)+smallRate)
        *(Foam::pow(krn,2)*dkrwdS/mun + Foam::pow(krw,2)*dkrndS/muw)
        /(muw/mun*Foam::pow(krn,2) + 2*krn*krw + mun/muw*Foam::pow(krw,2));

    tensorField CFLCoats
    (
        runTime.deltaT()*dfw*fvc::surfaceSum(mag(phi))/porosity
    );

    // Capillarity
    CFLCoats += 
        (runTime.deltaT()/porosity)
        *2*mag(pcModel->dpcdS())
        *fvc::surfaceSum(Kf*mesh.magSf()/mag(delta))
        *(krn*krw/(muw*krn+mun*krw));

    CFLCoats /= mesh.V();

    CFLUse = cmptMax(gMax(CFLCoats));
    Info<< "Coats Number max: " << CFLUse << endl;
    maxDeltaTFact = maxCo/(CFLUse + SMALL);
      
}
