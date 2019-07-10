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
    Foam::blackoilPhase

Description
    A phase model for Black oil models, with FVF, mu and rho updates.

SourceFiles
    blackoilPhase.C

\*---------------------------------------------------------------------------*/

#ifndef blackoilPhase_H
#define blackoilPhase_H

#include "phase.H"
#include "fvcFlux.H"
#include "FVFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class blackoilPhase Declaration
\*---------------------------------------------------------------------------*/

class blackoilPhase
:
public phase
{

protected:
    
    //- Phase density at standard conditions
    dimensionedScalar rhoSc_;

    //- Phase density field
    volScalarField rho_;

    //- Phase dynamic viscosity
    volScalarField mu_;

    //- Phase Formation Volume Factor
    volScalarField rFVF_;

    //- (drFVF|dP) field
    volScalarField drFVFdP_;

    //- FVF Model
    autoPtr<FVFModel> BModel_;

    //- Phase flux
    autoPtr<surfaceScalarField> phiPtr_;

public:

    // Constructors

        //- Construct from components
        blackoilPhase
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

    // Selectors

        //- Return a reference to incompressible phase
        static autoPtr<blackoilPhase> New
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~blackoilPhase();

    // Member Functions
    
        //- Return dynamic viscosity
        volScalarField& mu()
        {
            return mu_;
        }

        //- Return dynamic viscosity
        const volScalarField& mu() const
        {
            return mu_;
        }

        const dimensionedScalar standardRho() const
        {
            return rhoSc_;
        }

        const dimensionedScalar& standardRho()
        {
            return rhoSc_;
        }

        //- Return density
        volScalarField& rho()
        {
            return rho_;
        }

         //- Return density
        const volScalarField& rho() const
        {
            return rho_;
        }
        
        //- Return const rFVF
        volScalarField& rFVF() 
        {
            return rFVF_;
        }

        //- Return const drFVFdP
        volScalarField& drFVFdP() 
        {
            return drFVFdP_;
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

        //- Update FVF from Models
        void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
