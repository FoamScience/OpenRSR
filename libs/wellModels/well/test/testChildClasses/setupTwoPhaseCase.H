// The wellProperties dict
IOdictionary wellProperties
(
    IOobject
    (
        "wellProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
// The transportProperties dict (with rho and mu)
dictionary transportProperties;

// A sub dict with rho - mu as dimensionedScalars
dictionary waterDict;
waterDict.add("rho", dimensionedScalar(word("water"), dimMass/dimVolume, 1000));
waterDict.add("mu", dimensionedScalar(word("water"), dimMass/dimLength/dimTime, 1e-3));
// A sub dict with rho - mu as dimensionedScalars
dictionary oilDict;
oilDict.add("rho", dimensionedScalar(word("oil"), dimMass/dimVolume, 700));
oilDict.add("mu", dimensionedScalar(word("oil"), dimMass/dimLength/dimTime, 1e-5));

// Add the subdict to the parent transportDict
transportProperties.add(word("water"), waterDict);
// Add the subdict to the parent transportDict
transportProperties.add(word("oil"), waterDict);

// The pressure field
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("p", dimPressure, 0.0)
);

// Two phases to work with
incompressibleFluid oil(word("oil"), transportProperties, mesh);
incompressibleFluid water(word("water"), transportProperties, mesh);
