// This should belong to an Impes algo. class
scalar dSmax(runTime.controlDict().lookupOrDefault<scalar>("dSmax",0.));

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Collect phase names
wordList phaseNames(transportProperties.lookup("phases"));


// Phase objects
incompressibleFluid canPhase(phaseNames[0], transportProperties, mesh);
incompressibleFluid noncanPhase(phaseNames[1], transportProperties, mesh);

// Convenience references
volVectorField& Uc = canPhase.U();
surfaceScalarField& phic = canPhase.phi();
volVectorField& Un = noncanPhase.U();
surfaceScalarField& phin = noncanPhase.phi();

const uniformDimensionedScalarField& rhoc = canPhase.rho();
const uniformDimensionedScalarField& muc = canPhase.mu();
const uniformDimensionedScalarField& rhon = noncanPhase.rho();
const uniformDimensionedScalarField& mun = noncanPhase.mu();

// Relative permeability model
autoPtr<base2PhasesKrModel> krModel =
    base2PhasesKrModel::New("krModel", transportProperties, mesh);

// Capillary pressure model
autoPtr<base2PhasesPcModel> pcModel =
    base2PhasesPcModel::New("pcModel", transportProperties, mesh);

Info<< "Reading field porosity\n" << endl;
volScalarField porosity
(
    IOobject
    (
        "porosity",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field K\n" << endl;
volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Absolute Permeability interpolation to faces
surfaceScalarField Kf = fvc::interpolate(K,"K");

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Un + Uc
);

// Fix a bug with foam-extend-4.0
mesh.schemesDict().setFluxRequired(p.name());

// Hard-coded Diagonal solver for saturation
// I don't see the point of specifying this in fvSolution
dictionary SwSolver = dictionary();
SwSolver.add("solver","diagonal");
SwSolver.add("relTol",0.0);
Info << SwSolver << endl;
