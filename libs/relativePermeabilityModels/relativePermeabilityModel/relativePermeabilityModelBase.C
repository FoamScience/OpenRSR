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

#include "zeroGradientFvPatchField.H"

namespace Foam 
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int nPhases>
relativePermeabilityModelBase<nPhases>::relativePermeabilityModelBase
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
    :
    name_(name),
    transportProperties_(transportProperties),
    mesh_(mesh),
    kr_
    (
        IOobject
        (
            "krTable",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedKrVector<nPhases>( "kr", dimless, krVector<nPhases>::zero),
        zeroGradientFvPatchField<vector>::typeName
    ),
    dkrdS_
    (
        IOobject
        (
            "dkrdSTable",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedKrVector<nPhases>( "dkrdS", dimless, krVector<nPhases>::zero),
        zeroGradientFvPatchField<vector>::typeName
    )
{
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<int nPhases>
autoPtr<relativePermeabilityModelBase<nPhases> >
relativePermeabilityModelBase<nPhases>::New
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
{
    // Get the name from the dictionary.
    const word modelType(transportProperties.lookup("relativePermeabilityModel"));

    // Get the RTS Table via the global object.  
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);
    // If the constructor pointer is not found in the table.
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "relativePermeabilityModelBase::New(const dictionary&)"
        )   << "Unknown relativePermeabilityModelBase type "
            << modelType << nl << nl
            << "Valid relativePermeabilityModelBases are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr< relativePermeabilityModelBase<nPhases> > 
        (cstrIter()(name, transportProperties, mesh));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
