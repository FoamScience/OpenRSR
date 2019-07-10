/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "wellModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(wellModel, 0);
defineRunTimeSelectionTable(wellModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wellModel::wellModel
(
    const word& name,
    const dictionary& wellsProperties,
    const fvMesh& mesh,
    const HCmixture<blackoilPhase>& mixture,
    const relativePermeabilityModel& krModel
)
:
    name_(name),
    wellsProperties_(wellsProperties),
    mesh_(mesh),
    p_(mesh.lookupObject<volScalarField>("p")),
    mixture_(mixture),
    krModel_(krModel),
    wells_(well::wellsList),
    eoSource_
    (
        IOobject
        (
            "oilExpSource",
            p_.time().timeName(),
            p_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p_.mesh(),
        dimensionedScalar("",dimless/dimTime,0.0)
    ),
    ioSource_
    (
        IOobject
        (
            "oilImpSource",
            mesh.time().timeName()+"/"+typeName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("",dimless/dimTime/dimPressure,0.0)
    ),
    ewSource_
    (
        IOobject
        (
            "waterExpSource",
            mesh.time().timeName()+"/"+typeName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("",dimless/dimTime,0.0)
    ),
    iwSource_
    (
        IOobject
        (
            "waterImpSource",
            mesh.time().timeName()+"/"+typeName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("",dimless/dimTime/dimPressure,0.0)
    ),
    timeForDt_(VGREAT),
    WJw_(p_, dimless/dimTime)
{
    well::readWells(wellsProperties_, mesh_, mixture_);
    //forAll(wells_, welli)
    //{
    //    wells_[welli].correct();
    //    const labelList& cells = wells_[welli].wellSet().toc();

    //    forAll(cells, celli)
    //    {
    //        eoSource_[cells[celli]] = wells_[welli].eoSource()[celli];
    //        ewSource_[cells[celli]] = wells_[welli].ewSource()[celli];
    //        ioSource_[cells[celli]] = wells_[welli].ioSource()[celli];
    //        iwSource_[cells[celli]] = wells_[welli].iwSource()[celli];
    //    }
    //}
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::autoPtr<Foam::wellModel> Foam::wellModel::New
(
    const word& name,
    const dictionary& wellsProperties,
    const fvMesh& mesh,
    const HCmixture<blackoilPhase>& mixture,
    const relativePermeabilityModel& krModel
)
{
    const word wellManagement(wellsProperties.lookup("wellModel"));

    Info<< "Selecting wellModel : " << wellManagement << "\n" << endl;

    dictionaryConstructorTable::iterator cstrIter =
      dictionaryConstructorTablePtr_->find(wellManagement);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
      {
          FatalErrorIn("Foam::autoPtr<Foam::wellModel> Foam::wellModel::New\
                        (\
                            const word& name,\
                            const dictionary& wellsProperties,\
                            const fvMesh& mesh,\
                            const volScalarfield& p\
                        )" )
              << "Unknown wellManagement model "
              << wellManagement << nl << nl
              << "Valid wellManaegement models are : " << endl
              << dictionaryConstructorTablePtr_->sortedToc()
              << exit(FatalError);
      }

    return autoPtr<wellModel>
      (cstrIter()(name, wellsProperties, mesh, mixture, krModel));
}

// ************************************************************************* //
