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
    Foam::relativePermeabilityModels::krTabulated

Description
    Read evolution of the relative permeability and its derivative in respect 
    to the canonical phase.

SourceFiles
    krTabulated.C

\*---------------------------------------------------------------------------*/

#ifndef krTabulated_H
#define krTabulated_H

#include "phase.H"
#include "relativePermeabilityModelBase.H"
#include "relativePermeabilityModelBasesFwd.H"
#include "interpolationTable.H"
#include "vector4.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------*\
        Class krTabulated Declaration
\*---------------------------------------------------------------------------*/

class krTabulated
:
    public base2PhasesKrModel
{
    //- Private Data Members

        //- Coefficients dictionary
        dictionary coeffsDict_;

        //- Associated phase names
        wordList phaseNames_;

        //- Canonical phase name
        word canonicalPhase_;

        //- Const-ref to the canonical saturation
        const volScalarField& S_;

        //- Kr table
        interpolationTable<vector4> krSeries_;



public:

    //- Runtime type information
    TypeName("tabulated");

    // Constructors

        //- Construct from components
        krTabulated
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

        //- Destructor
        ~krTabulated() {}

    // Public Member Functions

        //- Return associated phase names as a list
        const wordList phaseNames() const 
        {
            return phaseNames_;
        }

        //- Return canonical phase name
        const word canonicalPhase() const
        {
            return canonicalPhase_;
        }

        //- Return phase index
        label phaseIndex(const word& phaseName) const
        {
            return (phaseName == canonicalPhase_ ? 0 : 1);
        }

        //- Correct the relative permeabilities
        void correct();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
