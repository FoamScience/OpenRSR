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
    Foam::blackoilWellModels::Peaceman

Description
    A class implementing a "workaround" for removing singularity when 
    solving for wells BHP; Taking the most simple case. 

SourceFiles
    Peaceman.C

\*---------------------------------------------------------------------------*/

#ifndef Peaceman_H
#define Peaceman_H

#include "wellModel.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blackoilWellModels
{

/*---------------------------------------------------------------------------*\
                  Class Peaceman Declaration
\*---------------------------------------------------------------------------*/

class Peaceman
:
    public wellModel
{

    // Private Data Members 

public:

    //- Runtime type information
    TypeName("Peaceman");

    // Constructors

        //- Construct from components
        Peaceman
        (
            const word& name,
            const dictionary& wellsProperties,
            const fvMesh& mesh,
            const HCmixture<blackoilPhase>& mixture,
            const relativePermeabilityModel& krModel
        );

    //- Destructor
    ~Peaceman() {}

    // Member Functions
    
        //- Correct well model
        void correct();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wellModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
