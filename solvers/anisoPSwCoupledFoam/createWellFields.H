
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
volScalarField& ewSource = wModel->ExplicitWaterSource();
volScalarField& eoSource = wModel->ExplicitOilSource();
volScalarField& iwSource = wModel->ImplicitWaterSource();
volScalarField& ioSource = wModel->ImplicitOilSource();
