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
    Foam::bubblePointModels::blackOilFlash

Description
    Attempting a simple BlackOil Flash operation to detect third phase formation
    ( Variable bubble point pressure ).

SourceFiles
    blackOilFlash.C

\*---------------------------------------------------------------------------*/

#ifndef blackOilFlash_H
#define blackOilFlash_H

#include "bubblePointModel.H"
#include "interpolationTable.H"
#include "vector2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace bubblePointModels
{

/*---------------------------------------------------------------------------*\
        Class blackOilFlash Declaration
\*---------------------------------------------------------------------------*/

class blackOilFlash
:
    public bubblePointModel
{
    //- Private Data Members
        
        //- Pb coeffs sub-dict
        dictionary PbCoeffs_;

        //- K-values fields
        volScalarField Ko_, Kg_;

        //- Solubility data series
        interpolationTable<vector2> RsSeries_;

        //- Standard molar densities for gas and oil
        dimensionedScalar rhoMost_,rhoMgst_;

        //- Oil & gas solubility data
        volScalarField Ro_, Rg_;

        //- Oil & gas molar solubility
        volScalarField ro_, rg_;

        //- Oil composition in liquid phase
        volScalarField xo_;

        //- Oil composition in gas phase
        volScalarField yo_;

        //- Gas composition in liquid phase
        volScalarField xg_;

        //- Gas composition in gas phase
        volScalarField yg_;

        
public:

    //- Runtime type information
    TypeName("blackOilFlash");

    // Constructors

    //- Construct from components
    blackOilFlash
    (
        const word& name,
        const dictionary& transportProperties,
        const volScalarField& p
    );

    //- Destructor
    ~blackOilFlash(){}

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
