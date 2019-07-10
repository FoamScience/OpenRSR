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

#include "capillaryPressureModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(capillaryPressureModel, 0);
defineRunTimeSelectionTable(capillaryPressureModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillaryPressureModel::capillaryPressureModel
(
    const word& name,
    const dictionary& transportProperties,
    const phase& wettingPhase,
    const phase& nonWettingPhase,
    const word& keyword
)
    :
    name_(name),
    transportProperties_(transportProperties),
    Sw_(wettingPhase.alpha()),
    pc_
    (
        IOobject
        (
            "pc",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sw_.mesh(),
        dimensionedScalar
        (
            "",
            dimensionSet(1,-1,-2,0,0,0,0),
            0.0
        )
    ),
    dpcdS_
    (
        IOobject
        (
            "dpcdS",
            Sw_.time().timeName(),
            Sw_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Sw_.mesh(),
        dimensionedScalar
        (
            "",
            dimensionSet(1,-1,-2,0,0,0,0),
            0.0
        )
    )
{}

Foam::autoPtr<Foam::capillaryPressureModel> Foam::capillaryPressureModel::New
(
    const word& name,
    const dictionary& transportProperties,
    const phase& wettingPhase,
    const phase& nonWettingPhase,
    const word& keyword
)
{
  const word modelType(transportProperties.lookup(keyword));

  Info<< "Selecting Capillary Pressure model:  " << modelType << "\n" << endl;

  dictionaryConstructorTable::iterator cstrIter =
    dictionaryConstructorTablePtr_->find(modelType);

  if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("Foam::autoPtr<Foam::capillaryPressureModel> Foam::capillaryPressureModel::New\
                    (\
                        const word& name,\
                        const dictionary& transportProperties,\
                        const phase& wettingPhase,\
                        const phase& nonWettingPhase,\
                        const word& keyword\
                     ) ")
	        << "Unknown capillaryPressureModel type "
	        << modelType << nl << nl
	        << "Valid capillaryPressureModels are : " << endl
	        << dictionaryConstructorTablePtr_->sortedToc()
	        << exit(FatalError);
    }

  return autoPtr<capillaryPressureModel>
    (cstrIter()(name, transportProperties, wettingPhase, nonWettingPhase));
}
// ************************************************************************* //
