// Create Time object
Foam::Time runTime
(
    Foam::Time::controlDictName,
    ".",
    "",
    "system",
    "constant",
    false
);

// Create Mesh object
Foam::fvMesh mesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
