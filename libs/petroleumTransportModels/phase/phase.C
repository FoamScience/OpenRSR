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

#include "phase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phase::phase
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
:  
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    phaseDict_(transportProperties.subDict(name)),
    mesh_(mesh),
    U_
    (
        IOobject
        (
            name+".U",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alpha_
    (
        IOobject
        (
            name+".alpha",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(name+".alpha", dimless, 1.0)
    ),
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
           linearInterpolate(U_) & mesh.Sf(),
           "fixedValue"
        )
    )
{
}
