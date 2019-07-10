
// Wells Dictionary
IOdictionary wellsProperties
(
    IOobject
    (
        "wellsProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Wells model
autoPtr<wellModel> wModel = wellModel::New("wModel", wellsProperties, mesh, mixture, krModel);
wModel->correct();

// Well Sources & Sinks
const volScalarField& ewSource = wModel->ExplicitWaterSource();
const volScalarField& eoSource = wModel->ExplicitOilSource();
const volScalarField& iwSource = wModel->ImplicitWaterSource();
const volScalarField& ioSource = wModel->ImplicitOilSource();

