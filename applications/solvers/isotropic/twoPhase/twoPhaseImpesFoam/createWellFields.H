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
