{

    // Plan to support water aquifer entries
    //
    //
    //// Store old values and update BCs
    //Un.correctBoundaryConditions();
    //Uw.correctBoundaryConditions();
    //// Set refs for phi if U is a fixedValue
    //forAll(mesh.boundary(),patchi)
    //{
    //    if (isA< fixedValueFvPatchField<vector> >(Un.boundaryField()[patchi]))
    //    {
    //        phin.boundaryField()[patchi] = Un.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
    //    }
    //    if (isA< fixedValueFvPatchField<vector> >(Uw.boundaryField()[patchi]))
    //    {
    //        phiw.boundaryField()[patchi] = Uw.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
    //    }
    //}


    volScalarField wellSource
    (
        IOobject
        (
            "wellSource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless/dimTime, 0),
        "zeroGradient"
    );

    forAll(wellSource, ci)
    {
        wellSource[ci] = wModel->explicitSource(canPhase.name())[ci];
    }

    fvScalarMatrix alphaEqn
    (
        porosity*fvm::ddt(canPhase.alpha())
        + fvc::div(phic)
    );
    alphaEqn.source() += wModel->explicitSource(canPhase.name());

    //// Use the hardcoded diagonal solver, disable annoying output
    //label oldDebugLevel = blockLduMatrix::debug();
    //blockLduMatrix::debug = 0;
    //alphaEqn.solve(SwSolver);
    //blockLduMatrix::debug = oldDebugLevel;
    ////lduMatrix::debug = oldDebugState;
    alphaEqn.solve();

    Info << canPhase.name() << " saturation: Min = "
        << gMin(canPhase.alpha()) << " Max = " << gMax(canPhase.alpha()) <<  endl;
}
