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

#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relativePermeabilityModel, 0);
    defineRunTimeSelectionTable(relativePermeabilityModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModel::relativePermeabilityModel
(
    const word& name,
    const dictionary& transportProperties,
    const phase& phase1,
    const phase& phase2,
    const word& keyword
)
    :
    name_(name),
    transportProperties_(transportProperties),
    phase1_(phase1),
    phase2_(phase2),
    S_(phase1.alpha()),
    kr1_
    (
        IOobject
        (
            phase1_.name()+".kr",
            S_.time().timeName(),
            S_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S_.mesh(),
        dimensionedScalar("kr", dimless, 0.0)
    ),
    kr2_
    (
        IOobject
        (
            phase2_.name()+".kr",
            S_.time().timeName(),
            S_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S_.mesh(),
        dimensionedScalar("kr", dimless, 0.0)
    ),
    dkr1dS_
    (
        IOobject
        (
            phase1_.name()+".dkrdS",
            S_.time().timeName(),
            S_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S_.mesh(),
        dimensionedScalar("dkrdS", dimless, 0.0)
    ),
    dkr2dS_
    (
        IOobject
        (
            phase2_.name()+".dkrdS",
            S_.time().timeName(),
            S_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S_.mesh(),
        dimensionedScalar("dkrdS", dimless, 0.0)
    )
{}

Foam::autoPtr<Foam::relativePermeabilityModel> Foam::relativePermeabilityModel::New
(
    const word& name,
    const dictionary& transportProperties,
    const phase& phase1,
    const phase& phase2,
    const word& keyword
)
{
    const word modelType(transportProperties.lookup(keyword));

    Info<< "Selecting relativePermeability model: " 
        << modelType << "\n" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "Foam::autoPtr<Foam::relativePermeabilityModel> Foam::relativePermeabilityModel::New\
            (\
                const word& name,\
                const dictionary& transportProperties,\
                const phase& phase1,\
                const phase& phase2\
            ) ")   
            << "Unknown relativePermeabilityModel type "
            << modelType << nl << nl
            << "Valid relativePermeabilityModels are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<relativePermeabilityModel>
        (cstrIter()(name, transportProperties, phase1, phase2));
}

// ************************************************************************* //
