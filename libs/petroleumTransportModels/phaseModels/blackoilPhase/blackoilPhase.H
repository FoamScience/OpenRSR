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
    A basic class for Black oil phases, with FVF, mu and rho updates.

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

    //- FVF Model
    autoPtr<FVFModel> BModel_;

public:

    // Constructors

        //- Construct from components
        blackoilPhase
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

        //- Destructor
        virtual ~blackoilPhase(){}

    // Public Member Functions

        //- Return a ref to dynamic viscosity
        volScalarField& mu()
        {
            return mu_;
        }

        //- Return a const-ref to dynamic viscosity
        const volScalarField& mu() const
        {
            return mu_;
        }

        //- Return a const-ref to standard density
        const dimensionedScalar& standardRho() const
        {
            return rhoSc_;
        }

        //- Return ref to density field
        volScalarField& rho()
        {
            return rho_;
        }

         //- Return const-ref to density field
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return const 1|FVF
        const volScalarField& rFVF() const
        {
            if (BModel_.valid())
            {
                return BModel_->rFVF();
            } else {
                FatalErrorIn("Trying to access rFVF field from a null pointer.") 
                    << exit(FatalError);
            }
        }

        //- Return const drFVFdP ref
        const volScalarField& drFVFdP()
        {
            if (BModel_.valid())
            {
                return BModel_->drFVFdP();
            } else {
                FatalErrorIn("Trying to access drFVF|dP field from a null pointer.") 
                    << exit(FatalError);
            }
        }

        //- Update FVF from Models
        void correct();

        //- Write data to ostream
        virtual bool writeData(Ostream& os) const 
        {
            notImplemented("bool writeData(Ostream& os) const");
            return 0;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
