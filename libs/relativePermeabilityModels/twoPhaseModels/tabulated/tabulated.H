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
    Foam::relativePermeabilityModels::tabulated

Description
    Tabular description for two-phase Kr as a function of a chosen saturation.

SourceFiles
    tabulated.C

\*---------------------------------------------------------------------------*/

#ifndef tabulated_H
#define tabulated_H

#include "relativePermeabilityModel.H"
#include "interpolationTable.H"
#include "vector4.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------*\
        Class tabulated Declaration
\*---------------------------------------------------------------------------*/

class tabulated
:
    public relativePermeabilityModel
{
    //- Private Data Members
        
        //- Relative permeability table
        interpolationTable<vector4> krSeries_;

        
public:

    //- Runtime type information
    TypeName("tabulated");

    // Constructors

    //- Construct from components
    tabulated
    (
        const word& name,
        const dictionary& transportProperties,
        const phase& phase1,
        const phase& phase2
    );

    //- Destructor
    ~tabulated()
        {}

    // Member Functions

    //- Correct the relative permeabilities
    void correct()
        {

            // Set internal field
            forAll(kr1_.internalField(), celli)
            {
                kr1_[celli] = krSeries_(S_.internalField()[celli]) [0];
                kr2_[celli] = krSeries_(S_.internalField()[celli]) [1];
                dkr1dS_[celli] = krSeries_(S_.internalField()[celli]) [2];
                dkr2dS_[celli] = krSeries_(S_.internalField()[celli]) [3];
            }

            forAll(kr1_.boundaryField(), patchi)
            {
                forAll(kr1_.boundaryField()[patchi], facei)
                {
                    scalar Satf =  S_.boundaryField()[patchi][facei];
                    kr1_.boundaryField()[patchi][facei] = krSeries_(Satf)[0];
                    kr2_.boundaryField()[patchi][facei] = krSeries_(Satf)[1];
                    dkr1dS_.boundaryField()[patchi][facei] = krSeries_(Satf)[2];
                    dkr2dS_.boundaryField()[patchi][facei] = krSeries_(Satf)[3];
                }
            }
            kr1_.correctBoundaryConditions();
            kr2_.correctBoundaryConditions();
            dkr1dS_.correctBoundaryConditions();
            dkr2dS_.correctBoundaryConditions();
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
