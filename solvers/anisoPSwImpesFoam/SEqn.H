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
  
    fvScalarMatrix SwEqn
        (
            porosity*fvm::ddt(phasew.alpha()) 
            + fvc::div(phiw) 
            // Wells contribution
            + ewSource+iwSource*p
        );

    // Use the hardcoded diagonal solver, disable annoying output
    label oldDebugLevel = blockLduMatrix::debug();
    blockLduMatrix::debug = 0;
    SwEqn.solve(SwSolver);
    blockLduMatrix::debug = oldDebugLevel;
    //lduMatrix::debug = oldDebugState;

    Info << "Saturation w " << " Min = " << gMin(phasew.alpha()) << " Max = " << gMax(phasew.alpha()) <<  endl;
}
