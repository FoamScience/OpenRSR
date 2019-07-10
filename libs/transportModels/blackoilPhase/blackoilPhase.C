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

Foam::blackoilPhase::blackoilPhase
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
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        rhoSc_,
        "zeroGradient"
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
    rFVF_(
        IOobject
        (
            name+".rFVF",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("",dimless, 0.0),
        "zeroGradient"
    ),
    drFVFdP_(
        IOobject
        (
            name+".drFVFdP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("",dimless/dimPressure, 0.0),
        "zeroGradient"
    ),
    BModel_
    ( 
        FVFModel::New(
           name_+".FVFModel", 
           phaseDict_,
           mesh_.objectRegistry::lookupObject<volScalarField>("p")
        )
    )
{
    //wordList phiBCTypes
    //(
    //    U_.boundaryField().size(),
    //    calculatedFvPatchScalarField::typeName // or just "calculated"
    //);
    
    phiPtr_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                name_+".phi",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            Foam::linearInterpolate(U_) & mesh.Sf(),
            calculatedFvPatchScalarField::typeName
        )
    );

    BModel_->correct();
    rFVF_ = BModel_->rFVF();
    drFVFdP_ = BModel_->drFVFdP();
}


Foam::autoPtr<Foam::blackoilPhase> Foam::blackoilPhase::New
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
{
    return autoPtr<blackoilPhase>
    (
        new blackoilPhase(name, transportProperties, mesh)
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blackoilPhase::~blackoilPhase()
{}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

void Foam::blackoilPhase::correct()
{
    // Run the model's correct() method to calculate FVF = f(P)
    BModel_->correct();
    rFVF_ = BModel_->rFVF();
    drFVFdP_ = BModel_->drFVFdP();

    rho_ = rhoSc_*rFVF_;
    rho_.correctBoundaryConditions();
}

// ************************************************************************* //
