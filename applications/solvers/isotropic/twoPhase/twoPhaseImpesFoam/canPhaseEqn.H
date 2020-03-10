{
    // Plan to support water aquifer entries
    //
    //
    //// Store old values and update BCs
    //Uc.correctBoundaryConditions();
    //Un.correctBoundaryConditions();
    //// Set refs for phi if U is a fixedValue
    //forAll(mesh.boundary(),patchi)
    //{
    //    if (isA< fixedValueFvPatchField<vector> >(Un.boundaryField()[patchi]))
    //    {
    //        phin.boundaryField()[patchi] = Un.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
    //    }
    //    if (isA< fixedValueFvPatchField<vector> >(Uc.boundaryField()[patchi]))
    //    {
    //        phic.boundaryField()[patchi] = Uc.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
    //    }
    //}

    fvScalarMatrix canPhaseEqn
    (
        porosity*fvm::ddt(canPhase.alpha())
        + fvc::div(phic)
    );
    // Explicit well contribution
    canPhaseEqn.source() += wModel->explicitSource(canPhase.name());

    // Use the hardcoded diagonal solver, disable annoying output
    if (!mesh.solutionDict().subDict("solvers").found(canPhase.name()+".alpha"))
    {
        label oldDebugLevel = blockLduMatrix::debug();
        blockLduMatrix::debug = 0;
        canPhaseEqn.solve(SwSolver);
        blockLduMatrix::debug = oldDebugLevel;
        if(fvSchemes::debug)
        {
            Info<< "Hard-coded scheme for the saturation equation:"
                << nl << SwSolver << endl;
        }
    } else {
        canPhaseEqn.solve();
    }

    Info << canPhase.name() << " saturation: Min = "
         << gMin(canPhase.alpha()) << " Max = "
         << gMax(canPhase.alpha()) <<  endl;
}
