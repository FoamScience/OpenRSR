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
    name_(name),
    phaseDict_(transportProperties.subDict(name)),
    mesh_(mesh),
    U_
    (
        IOobject
        (
            name+".U",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    alpha_
    (
        IOobject
        (
            name+".alpha",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(name+".alpha", dimless, 1.0)
    )
{
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

//Foam::autoPtr<Foam::phase> Foam::phase::New
//(
//    const word& name,
//    const dictionary& transportProperties,
//    const fvMesh& mesh
//)
//{
//    return autoPtr<phase>
//    (
//        new phase(name, transportProperties, mesh)
//    );
//}
