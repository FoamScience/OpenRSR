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
capillaryPressureModelBase<nPhases>::capillaryPressureModelBase
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
    :
    name_(name),
    transportProperties_(transportProperties),
    mesh_(mesh),
    pc_
    (
        IOobject
        (
            "pcTable",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<cmptT<nPhases> >
            ( "pc", dimless, pcTableReturn<nPhases>::zero()),
        zeroGradientFvPatchField<vector>::typeName
    ),
    dpcdS_
    (
        IOobject
        (
            "dpcdSTable",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<cmptT<nPhases> >
            ( "dpcdS", dimless, pcTableReturn<nPhases>::zero()),
        zeroGradientFvPatchField<vector>::typeName
    )
{
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<int nPhases>
autoPtr<capillaryPressureModelBase<nPhases> >
capillaryPressureModelBase<nPhases>::New
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
{
    // Get the name from the dictionary.
    const word modelType(transportProperties.lookup("capillaryPressureModel"));

    // Get the RTS Table via the global object.  
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);
    // If the constructor pointer is not found in the table.
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(__PRETTY_FUNCTION__)
            << "Unknown capillaryPressureModelBase type "
            << modelType << nl << nl
            << "Valid capillaryPressureModelBases are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr< capillaryPressureModelBase<nPhases> > 
        (cstrIter()(name, transportProperties, mesh));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
