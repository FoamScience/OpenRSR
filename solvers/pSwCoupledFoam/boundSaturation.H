{
    // Bound saturation field
    dimensionedScalar S1 = min(phasew.alpha());
    dimensionedScalar S2 = max(phasew.alpha());

    if (S1 < SMin || S2 > SMax)
    {
        Info<< "S: " << S1.value() << " " << S2.value()
            << ".  Bounding." << endl;

        phasew.alpha().max(SMin);
        phasew.alpha().min(SMax);
        // The solver should update BCs
        Info << gMax(krw) << tab << gMin(krw) << endl;
        Info << gMax(krn) << tab << gMin(krn) << endl;
    }
}
