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
    Foam::twoPhaseWells::verticalIsoWell

Description
    A vertical well specialization for isotropic solvers.

SourceFiles
    verticalIsoWell.C

\*---------------------------------------------------------------------------*/

#ifndef verticalIsoWell_H
#define verticalIsoWell_H

#include "wellBaseAndModel.H"
#include "wellBasesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace twoPhaseWells
{

/*---------------------------------------------------------------------------*\
                  Class verticalIsoWell Declaration
\*---------------------------------------------------------------------------*/

class verticalIsoWell
: public base2IsoWell
{

public:

    //- Runtime type information
    TypeName("verticalIsoWell");

    // Constructors

        //- Construct from components
        verticalIsoWell
        (
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh,
            const base2IsoWellModel& corrector
        );

    //- Destructor
    ~verticalIsoWell() {}

    // Member Functions

        //- Write Well data to Ostream
        virtual bool writeData(Ostream&) const;

        //- Pre-Setup for vertical well correction
        virtual void preCorrect();

        //- Post-Cleanup on vertical well correction
        virtual void postCorrect();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}; // End namespace twoPhaseWells

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
