//  Create fields necessary to drive the simulation

//- Max change in canonical phase saturation through a timeStep
scalar dSmax(runTime.controlDict().lookupOrDefault<scalar>("dSmax", 0.01));

Info<< nl << "Reading field p\n" << endl;
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

// Get Phase names from transportProperties
wordList phaseNames;
transportProperties.lookup("phases") >> phaseNames;
if (phaseNames.size() != 2)
{
    FatalErrorIn("createFields.H")
        << "Expected 2 phase names in transportProperties.phases but"
        << phaseNames.size() << " incountered" << nl << exit(FatalError);
}

//- Two incompressible phase objects
phaseModels::incompressibleFluid canPhase(phaseNames[0], transportProperties, mesh);
phaseModels::incompressibleFluid noncanPhase(phaseNames[1], transportProperties, mesh);

// Convenience references for canonical (c) and noncanonical (n) phase parameters
volVectorField& Un = noncanPhase.U();
volVectorField& Uc = canPhase.U();
surfaceScalarField& phin = noncanPhase.phi();
surfaceScalarField& phic = canPhase.phi();

const uniformDimensionedScalarField& rhon = noncanPhase.rho();
const uniformDimensionedScalarField& rhoc = canPhase.rho();
const uniformDimensionedScalarField& mun = noncanPhase.mu();
const uniformDimensionedScalarField& muc = canPhase.mu();

// Relative permeability model
autoPtr<base2PhasesKrModel> krModel =
    base2PhasesKrModel::New("krModel", transportProperties, mesh);

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

// Interpolate Absolute Permeability to mesh faces
Info<< "Interpolating field K\n" << endl;
surfaceScalarField Kf = fvc::interpolate(K, "K");

// Sum of phase velocities in a cell
// Used to handle total flux!
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
    Uc + Un
);

// Require flux for p field here
// because foam-extend-4.0 won't read it from fvSchemes for some reason
mesh.schemesDict().setFluxRequired(p.name());
Info << nl;

dictionary SwSolver = dictionary();
SwSolver.add("solver","diagonal");
SwSolver.add("relTol",0.0);
