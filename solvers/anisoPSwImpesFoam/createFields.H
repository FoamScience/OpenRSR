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

// Blackoil phase objects within a mixtue
autoPtr<HCmixture<blackoilPhase> > mixture = HCmixture<blackoilPhase>::New("twoPhase",transportProperties,mesh);

// Access to phases by wettebility status, deprecated
blackoilPhase& phasen = mixture->phase("nonWetting");
blackoilPhase& phasew = mixture->phase("wetting");

// NonWetting phase convenient parameters
volVectorField& Un = phasen.U();
surfaceScalarField& phin = phasen.phi();
volScalarField& rhon = phasen.rho();
volScalarField& mun = phasen.mu();
// Wetting phase convenient parameters
volVectorField& Uw = phasew.U();
surfaceScalarField& phiw = phasew.phi();
volScalarField& rhow = phasew.rho();
volScalarField& muw = phasew.mu();

// Relative permeability model 
autoPtr<relativePermeabilityModel> krModel = 
    relativePermeabilityModel::New("krModel", transportProperties, phasew, phasen);

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
volTensorField K
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
surfaceTensorField Kf = fvc::interpolate(K,"K");

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
    Un + Uw
);

// Capillary pressure
autoPtr<capillaryPressureModel> pcModel = 
    capillaryPressureModel::New("pcModel", transportProperties, phasew, phasen);

// Fix a bug with foam-extend-4.0
mesh.schemesDict().setFluxRequired(p.name());

// Hard-coded Diagonal solver for saturation
// I don't see the point of specifying this in fvSolution
dictionary SwSolver = dictionary();
SwSolver.add("solver","diagonal");
SwSolver.add("relTol",0.0);
Info << SwSolver << endl;
