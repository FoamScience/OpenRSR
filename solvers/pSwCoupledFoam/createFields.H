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
volScalarField p0 ("p0", p);

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
volScalarField& rFVFn = phasen.rFVF();
// Wetting phase convenient parameters
volVectorField& Uw = phasew.U();
surfaceScalarField& phiw = phasew.phi();
volScalarField& rhow = phasew.rho();
volScalarField& muw = phasew.mu();
volScalarField& rFVFw = phasew.rFVF();

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
volScalarField porosity_0 ("porosity_0", porosity);

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

Info<< "Reading field compressibility\n" << endl;
volScalarField compressibility
(
    IOobject
    (
        "compressibility",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("compressibility", dimless/dimPressure,
        transportProperties.found("rockCoeffs")
        ?  transportProperties.subDict("rockCoeffs").lookupOrDefault<scalar>("compressibility", 0.0)
        : 0.0
    ),
    "zeroGradient"
);

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
//dictionary SwSolver = dictionary();
//SwSolver.add("solver","diagonal");
//SwSolver.add("relTol",0.0);

// Coupled working variable
volVector2Field Swp
(
    IOobject
    (
        "Swp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector2("zero",dimless,vector2::zero)
);

