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
    Foam::FVFModels

Description
    A namespace for Formation Volume Factor models.

Class
    Foam::FVFModel

Description
    A generic abstract class for Formation Volume Factor.

SourceFiles
    FVFModel.C

\*---------------------------------------------------------------------------*/

#ifndef FVFModel_H
#define FVFModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class FVFModel Declaration
\*---------------------------------------------------------------------------*/

class FVFModel
{

protected:

    // Protected data

        //- Model Name 
        word name_;

        //- Hosting dict
        dictionary phaseDict_;

        //- Const ref to phase pressure field
        const volScalarField& p_; 

        //- 1|FVF field
        volScalarField rFVF_;

        //- d(1|FVF)/dP field
        volScalarField drFVFdP_;

        //- Disallow copy construct
        FVFModel(const FVFModel&);

        //- Disallow default bitwise assignment
        void operator=(const FVFModel&);


public:

    //- Runtime type information
    TypeName("FVFModel");


    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        FVFModel,
        dictionary,
        (
            const word& name,
            const dictionary& phaseDict,
            const volScalarField& p
        ),
        (name, phaseDict, p)
    );


    // Selectors

        //- Return a reference to the selected FVFModel model
        static autoPtr<FVFModel> New
        (
            const word& name,
            const dictionary& phaseDict,
            const volScalarField& p
        );


    // Constructors

        //- Construct from components
        FVFModel
        (
            const word& name,
            const dictionary& phaseDict,
            const volScalarField& p
        );


    //- Destructor
    virtual ~FVFModel()
        {}


    // Member Functions

        //- Return model hosting dictionary
        const dictionary& phaseDict() const
        {
            return phaseDict_;
        }

        //- Return 1|FVF
        const volScalarField& rFVF() const
        {
            return rFVF_;
        }

        //- Return d(1|FVF)/dP
        const volScalarField& drFVFdP() const
        {
            return drFVFdP_;
        }
          
        //- Correct the FVF
        virtual void correct() = 0; 

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
