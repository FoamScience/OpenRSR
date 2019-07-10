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
    Incompressible phase model. No native incompressible viscosity models
    support yet.

SourceFiles
    incompressibleFluid.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleFluid_H
#define incompressibleFluid_H

#include "phase.H"
#include "fvcFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
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
    dimensionedScalar rho_;

    //- Phase flux
    autoPtr<surfaceScalarField> phiPtr_;

    //- Phase dynamic viscosity
    dimensionedScalar mu_;

public:

    // Constructors

        //- Construct from components
        incompressibleFluid
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

    // Selectors

        //- Return a reference to incompressible phase
        static autoPtr<incompressibleFluid> New
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~incompressibleFluid();

    // Member Functions
    
        //- Return dynamic viscosity
        const dimensionedScalar& mu() const
        {
            return mu_;
        }

        //- Return density
        const dimensionedScalar& rho() const
        {
            return rho_;
        }
        
        //- Return const ref to flux
        const surfaceScalarField& phi() const
        {
            return phiPtr_();
        }
        
        //- Return ref to flux
        surfaceScalarField& phi()
        {
            return phiPtr_();
        }

        //- No need to correct stuff in incompressible flow
        void correct(){}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
