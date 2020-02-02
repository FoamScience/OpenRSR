/*---------------------------------------------------------------------------*\
  =========                 |
  am
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

#include "FVFModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FVFModel, 0);
    defineRunTimeSelectionTable(FVFModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FVFModel::FVFModel
(
    const word& name,
    const dictionary& phaseDict,
    const volScalarField& p
)
    :
    name_(name),
    phaseDict_(phaseDict),
    p_(p),
    rFVF_
    (
        IOobject
        (
            phaseDict.name()+".rFVF",
            p.mesh().time().timeName(),
            p.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        p.mesh(),
        dimensionedScalar("rFVF", dimless, 1.0),
        zeroGradientFvPatchField<vector>::typeName
    ),
    drFVFdP_
    (
        IOobject
        (
            phaseDict_.name()+".drFVFdP",
            p_.mesh().time().timeName(),
            p_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        p_.mesh(),
        dimensionedScalar("drFVFdP", dimless/dimPressure, 0.0),
        zeroGradientFvPatchField<vector>::typeName
    )
{}

Foam::autoPtr<Foam::FVFModel> Foam::FVFModel::New
(
    const word& name,
    const dictionary& phaseDict,
    const volScalarField& p
)
{
    const word modelType(phaseDict.lookup("FVFModel"));

    Info<< "Selecting FVF model : " << modelType
        << "\n" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(__PRETTY_FUNCTION__)
               << "Unknown FVFModel type "
                << modelType << nl << nl
                << "Valid FVFModels are : " << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
    }

    return autoPtr<FVFModel> (cstrIter()(name, phaseDict, p));
}

// ************************************************************************* //
