    // Read field bounds
    dictionary fieldBounds = mesh.solutionDict().subDict("fieldBounds");

    // Pressure bounds
    dimensionedScalar SMin("SMin", phasew.alpha().dimensions(), 0);
    dimensionedScalar SMax("SMax", phasew.alpha().dimensions(), 0);


    fieldBounds.lookup(phasew.alpha().name()) >> SMin.value() >> SMax.value();
