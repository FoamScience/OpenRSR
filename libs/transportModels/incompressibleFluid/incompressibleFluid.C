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

#include "incompressibleFluid.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleFluid::incompressibleFluid
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
:
    phase(name,transportProperties,mesh),
    rho_(phaseDict_.lookup("rho")),
    phiPtr_
    ( 
         new surfaceScalarField
         (
             IOobject
             (
                 name+".phi",
                 mesh.time().timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
             ),
             linearInterpolate(U_) & mesh.Sf()
         )
    ),
    mu_(phaseDict_.lookup("mu"))
{}


Foam::autoPtr<Foam::incompressibleFluid> Foam::incompressibleFluid::New
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
{
    return autoPtr<incompressibleFluid>
    (
        new incompressibleFluid(name, transportProperties, mesh)
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::incompressibleFluid::~incompressibleFluid()
{}

// ************************************************************************* //
