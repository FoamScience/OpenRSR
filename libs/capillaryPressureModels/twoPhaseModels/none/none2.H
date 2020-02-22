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

Class
    Foam::capillaryPressureModels::none

Description
    Disables capillary pressure modeling

SourceFiles
    none.C

\*---------------------------------------------------------------------------*/

#ifndef none_H
#define none_H

#include "phase.H"
#include "capillaryPressureModelBase.H"
#include "capillaryPressureModelBasesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace capillaryPressureModels
{

/*---------------------------------------------------------------------------*\
        Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public base2PhasesPcModel
{

public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        none
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

        //- Destructor
        ~none() {}

    // Public Member Functions

        //- Return associated phase names as a list
        const wordList phaseNames() const 
        {
            return wordList();
        }

        //- Return canonical phase name
        const word canonicalPhase() const
        {
            return word();
        }

        //- Correct the relative permeabilities
        virtual void correct() {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillaryPressureModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
