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
    Foam::incompressibleFluid

Description
    Incompressible phase model. No need for native incompressible viscosity models
    support (No meant to use in Pipeline calculations :) ).

SourceFiles
    incompressibleFluid.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleFluid_H
#define incompressibleFluid_H

#include "phase.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace phaseModels
{

/*---------------------------------------------------------------------------*\
               Class incompressibleFluid Declaration
\*---------------------------------------------------------------------------*/

class incompressibleFluid
:
public phase
{

protected:

    //- Phase density
    uniformDimensionedScalarField rho_;

    //- Phase dynamic viscosity
    uniformDimensionedScalarField mu_;

public:

    // Constructors

        //- Construct from components
        incompressibleFluid
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

        //- Destructor
        virtual ~incompressibleFluid(){}

    // Member Functions
    
        //- Return dynamic viscosity
        const uniformDimensionedScalarField& mu() const
        {
            return mu_;
        }

        //- Return density
        const uniformDimensionedScalarField& rho() const
        {
            return rho_;
        }
        
        //- Correct the velocity
        virtual void correct();

        //- Write data to ostream
        virtual bool writeData(Ostream& os) const 
        {
            notImplemented("bool writeData(Ostream& os) const");
            return 0;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
