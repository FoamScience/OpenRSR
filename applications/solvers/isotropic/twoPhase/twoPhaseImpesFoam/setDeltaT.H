{

    // Timestepping Saturation equation
    //volScalarField wellSource("wellSource", wModel->explicitSource(canPhase.name()));
    scalarField dSdt= mag(-fvc::div(phic)/*+wellSource+()*/)().internalField()
                      / porosity.internalField();
    scalar dtForS = dSmax/(gMax(dSdt)+SMALL);
    
    scalar CFLUse = -1;
    scalar maxDeltaTFact = -1;
    #include "impesCoatsNo.H"
    
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
    
    // Set the timeStep
    runTime.setDeltaT
    (
        min(
            dtForS,
            min(
                deltaTFact*runTime.deltaTValue(),
                maxDeltaT
            )
        )
    );
    
    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

// ************************************************************************* //
