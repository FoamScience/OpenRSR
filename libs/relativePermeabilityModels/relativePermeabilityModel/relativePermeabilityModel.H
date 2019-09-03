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

Namespace
    Foam::relativePermeabilityModels

Description
    A namespace for relative permeability models.

Class
    Foam::relativePermeabilityModel

Description
    A generic abstract class for tow-phase relative permeability models.

SourceFiles
    relativePermeabilityModel.C

\*---------------------------------------------------------------------------*/

#ifndef relativePermeabilityModel_H
#define relativePermeabilityModel_H

#include "volFields.H"
#include "HCmixture.H"
#include "phase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class relativePermeabilityModel Declaration
\*---------------------------------------------------------------------------*/

class relativePermeabilityModel
{


protected:

    // Protected data

        //- Model Name 
        word name_;

        //- Hosting dict
        dictionary transportProperties_;

        //- Refs to phases
        const phase& phase1_;
        const phase& phase2_;

        //- Const ref to reference saturation
        //  (Should be the incresaing one)
        const volScalarField& S_;

        //- Phase Relative permeabilities
        volScalarField kr1_,kr2_;

        //- Phase Relative permeability derivatives
        volScalarField dkr1dS_,dkr2dS_;

        //- Disallow copy construct
        relativePermeabilityModel(const relativePermeabilityModel&);

        //- Disallow default bitwise assignment
        void operator=(const relativePermeabilityModel&);


public:

    //- Runtime type information
    TypeName("relativePermeabilityModel");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        relativePermeabilityModel,
        dictionary,
        (
            const word& name,
            const dictionary& transportProperties,
            const phase& phase1,
            const phase& phase2
        ),
        (name, transportProperties, phase1, phase2)
    );


    // Selectors

        //- Return a reference to the selected relativePermeability model
        static autoPtr<relativePermeabilityModel> New
        (
            const word& name,
            const dictionary& transportProperties,
            const phase& phase1,
            const phase& phase2,
            const word& keyword = "relativePermeabilityModel"
        );


    // Constructors

        //- Construct from components
        relativePermeabilityModel
        (
            const word& name,
            const dictionary& transportProperties,
            const phase& phase1,
            const phase& phase2,
            const word& keyword = "relativePermeabilityModel"
        );


    //- Destructor
    virtual ~relativePermeabilityModel()
        {}


    // Member Functions

        //- Return model hosting dictionary
        const dictionary& transportProperties() const
        {
            return transportProperties_;
        }

        //- Return nonwetting relative permeability
        volScalarField& phase2Kr()
        {
            return kr2_;
        }

        //- Return wetting relative permeability
        volScalarField& phase1Kr()
        {
            return kr1_;
        }

        //- return derivative of nonwetting relative permeability
        volScalarField& phase2dKrdS()
        {
            return dkr2dS_;
        }   

        //- return derivative of nonwetting relative permeability
        volScalarField& phase1dKrdS()
        {
            return dkr1dS_;
        }

        //- Correct the relative permeabilities
        virtual void correct() = 0; 

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
