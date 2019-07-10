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
    Foam::FVFModels::tabulatedFVF

Description
    Tabular description of the relationship 1|FVF = f(P)

SourceFiles
    tabulatedFVF.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedFVF_H
#define tabulatedFVF_H

#include "FVFModel.H"
#include "vector2.H"
#include "interpolationTable.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FVFModels
{

/*---------------------------------------------------------------------------*\
    Class tabulatedFVF Declaration
\*---------------------------------------------------------------------------*/

class tabulatedFVF
:
    public FVFModel
{
    //- Private Data Members
        
        //- Table of pressure vs (rFVF drFVFdP)
        interpolationTable<vector2> rFVFseries_;

        
public:

    //- Runtime type information
    TypeName("tabulated");

    // Constructors

    //- Construct from components
    tabulatedFVF
    (
        const word& name,
        const dictionary& phaseDict,
        const volScalarField& p
    );

    //- Destructor
    ~tabulatedFVF()
        {}

    // Member Functions

    //- Correct the FVF
    void correct()
        {
            forAll(p_.internalField(), celli)
            {
                rFVF_[celli] = rFVFseries_(p_.internalField()[celli])[0];
                drFVFdP_[celli] = rFVFseries_(p_.internalField()[celli])[1];
            }
        
            // Make sure BC types are zeroGradient
            rFVF_.correctBoundaryConditions();
            drFVFdP_.correctBoundaryConditions();
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FVFModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
