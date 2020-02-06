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

#include "blackoilPhase.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModels::blackoilPhase::blackoilPhase
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
:
    phase(name, transportProperties, mesh),
    rhoSc_(phaseDict_.lookup("rhoSc")),
    rho_(
        IOobject
        (
            name+".rho",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        rhoSc_,
        zeroGradientFvPatchField<vector>::typeName
    ),
    mu_(
        IOobject
        (
            name+".mu",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh
    ),
    BModel_
    ( 
        FVFModel::New(
           name+".FVFModel", 
           transportProperties.subDict(name),
           mesh.objectRegistry::lookupObject<volScalarField>("p")
        )
    )
{
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

void Foam::phaseModels::blackoilPhase::correct()
{
    // Run the model's correct() method to calculate FVF = f(P)
    BModel_->correct();

    rho_ = rhoSc_*BModel_->rFVF();
    rho_.correctBoundaryConditions();
}

// ************************************************************************* //
