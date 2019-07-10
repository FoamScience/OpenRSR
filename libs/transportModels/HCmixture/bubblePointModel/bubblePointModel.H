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
    Foam::bubblePointModels

Description
    A namespace for bubble points detection classes

Class
    Foam::bubblePointModel

Description
    A generic abstract class for detection of bubble point pressure
    using BlackOil flash methods. 

SourceFiles
    bubblePointModel.C

\*---------------------------------------------------------------------------*/

#ifndef bubblePointModel_H
#define bubblePointModel_H

#include "volFields.H"
#include "phase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class bubblePointModel Declaration
\*---------------------------------------------------------------------------*/

class bubblePointModel
{

protected:

    // Protected data

        //- Model Name 
        word name_;

        //- Hosting dict
        dictionary transportProperties_;

        //- Const ref to pressure field
        const volScalarField& p_;

        //- Phases State:
        //  0 or 1: saturated
        //      -1: undersaturated
        volScalarField state_;

        //- Disallow copy construct
        bubblePointModel(const bubblePointModel&);

        //- Disallow default bitwise assignment
        void operator=(const bubblePointModel&);


public:

    //- Runtime type information
    TypeName("bubblePointModel");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        bubblePointModel,
        dictionary,
        (
            const word& name,
            const dictionary& transportProperties,
            const volScalarField& p
        ),
        (name, transportProperties, p)
    );


    // Selectors

        //- Return a reference to the selected relativePermeability model
        static autoPtr<bubblePointModel> New
        (
            const word& name,
            const dictionary& transportProperties,
            const volScalarField& p
        );


    // Constructors

        //- Construct from components
        bubblePointModel
        (
            const word& name,
            const dictionary& transportProperties,
            const volScalarField& p
        );


    //- Destructor
    virtual ~bubblePointModel()
        {}


    // Member Functions

        //- Return model hosting dictionary
        const dictionary& transportProperties() const
        {
            return transportProperties_;
        }

        //- Return phase state
        const volScalarField& state() const
        {
            return state_;
        }

          
        //- Correct the bubble point model
        virtual void correct() = 0; 

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
