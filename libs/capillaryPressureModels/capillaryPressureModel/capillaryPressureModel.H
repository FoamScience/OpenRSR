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
    Foam::capillaryPressureModels

Description
    A namespace for capillary pressure models.

Class
    Foam::capillaryPressureModel

Description
    An abstract class for capillary pressure models.

SourceFiles
    capillaryPressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef capillaryPressureModel_H
#define capillaryPressureModel_H

#include "volFields.H"
#include "phase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class capillaryPressureModel Declaration
\*---------------------------------------------------------------------------*/

class capillaryPressureModel
{

protected:

    // Protected data
    
        //- Model name
        word name_;

        //- Hosting dict
        dictionary transportProperties_;

        //- Const ref to wetting saturation
        const volScalarField& Sw_;

        //- Capillary pressure and its derivative 
        volScalarField pc_,dpcdS_;

    // Protected Member Functions

        //- Disallow copy construct
        capillaryPressureModel(const capillaryPressureModel&);

        //- Disallow default bitwise assignment
        void operator=(const capillaryPressureModel&);

public:

    //- Runtime type information
    TypeName("capillaryPressureModel");


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        capillaryPressureModel,
        dictionary,
        (
            const word& name,
            const dictionary& transportProperties,
            const phase& wettingPhase,
            const phase& nonWettingPhase
        ),
        (name, transportProperties, wettingPhase, nonWettingPhase)
    );


    // Selectors

        //- Return a reference to the selected capillary pressure model
        static autoPtr<capillaryPressureModel> New
        (
            const word& name,
            const dictionary& transportProperties,
            const phase& wettingPhase,
            const phase& nonWettingPhase,
            const word& keyword = "capillaryPressureModel"
        );


    // Constructors

        //- Construct from components
        capillaryPressureModel
        (
            const word& name,
            const dictionary& transportProperties,
            const phase& wettingPhase,
            const phase& nonWettingPhase,
            const word& keyword = "capillaryPressureModel"
        );


    //- Destructor
    virtual ~capillaryPressureModel() {}


    // Member Functions

        //- Return hosting dictionary
        const dictionary& transportProperties() const
        {
            return transportProperties_;
        }
          
        //- Return capillary pressure field
        tmp<volScalarField> pc() const
        {
            return pc_;
        }

        //- Return the derivative of capillary pressure
        tmp<volScalarField> dpcdS() const
        {
            return dpcdS_;
        }

        //- Correct the capillary pressure
        virtual void correct() = 0; 

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
