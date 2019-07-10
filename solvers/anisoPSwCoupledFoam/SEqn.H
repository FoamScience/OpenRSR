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
            fvm::ddt(wStorage, phasew.alpha()) 
            + fvc::div(phiwG) + fvc::div(phiPc)
            // Wells contribution
            + ewSource
        );

    fvScalarMatrix pInSw
        (
            fvm::laplacian(-Mwf,p)
            + fvm::Sp(iwSource,p)
        );

    resEqn.insertEquation(0,SwEqn);
    resEqn.insertEquationCoupling(0,1,pInSw);

}
