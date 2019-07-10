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
    Foam::bubblePointModels::constantPb

Description
    Fail the solver if fluid becomes saturated (constant bubble pressure)

SourceFiles
    constantPb.C

\*---------------------------------------------------------------------------*/

#ifndef constantPb_H
#define constantPb_H

#include "bubblePointModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace bubblePointModels
{

/*---------------------------------------------------------------------------*\
        Class constantPb Declaration
\*---------------------------------------------------------------------------*/

class constantPb
:
    public bubblePointModel
{
    //- Private Data Members
        
        //- Pb coeffs sub-dict
        dictionary PbCoeffs_;

        //- Pb to watch for
        dimensionedScalar Pb_;

        
public:

    //- Runtime type information
    TypeName("constantPb");

    // Constructors

    //- Construct from components
    constantPb
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& p
    );

    //- Destructor
    ~constantPb(){}

    // Member Functions

    //- Correct the phase state
    void correct();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace bubblePointModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
