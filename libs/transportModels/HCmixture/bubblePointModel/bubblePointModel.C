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

#include "bubblePointModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(bubblePointModel, 0);
defineRunTimeSelectionTable(bubblePointModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubblePointModel::bubblePointModel
(
    const word& name,
    const dictionary& transportProperties,
    const volScalarField& p
)
    :
    name_(name),
    transportProperties_(transportProperties),
    p_(p),
    state_
    (
        IOobject
        (
            "SaturationState",
            p_.time().timeName(),
            p_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        p_.mesh(),
        dimensionedScalar("", dimless, 0.0)
    )
{}

Foam::autoPtr<Foam::bubblePointModel> Foam::bubblePointModel::New
(
    const word& name,
    const dictionary& transportProperties,
    const volScalarField& p
)
{
    const word modelType(transportProperties.lookup("bubblePointModel"));

    Info<< "Selecting bubblePoint model: " 
        << modelType << "\n" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "Foam::autoPtr<Foam::bubblePointModel> Foam::bubblePointModel::New\
            (\
                const word& name,\
                const dictionary& transportProperties,\
                const volScalarField& p\
            ) ")   
            << "Unknown bubblePointModel type "
            << modelType << nl << nl
            << "Valid bubblePointModels are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<bubblePointModel>
        (cstrIter()(name, transportProperties, p));
}
// ************************************************************************* //
