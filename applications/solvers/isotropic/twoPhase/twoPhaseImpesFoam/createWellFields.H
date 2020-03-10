// Wells Dictionary
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

// Wells model
autoPtr<base2IsoWellModel> wModel =
	base2IsoWellModel::New("wModel", wellProperties, mesh);
wModel->correct();

// Well Sources & Sinks
//const volScalarField& ewSource = wModel->ExplicitWaterSource();
//const volScalarField& eoSource = wModel->ExplicitOilSource();
//const volScalarField& iwSource = wModel->ImplicitWaterSource();
//const volScalarField& ioSource = wModel->ImplicitOilSource();

