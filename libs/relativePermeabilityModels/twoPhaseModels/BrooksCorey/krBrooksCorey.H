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
    Foam::relativePermeabilityModels::krBrooksCorey

Description
    Modified two-phases Brooks-Corey relative permeability model:

    \rst

    .. math::

        k_{rw} = k_{rw, max}(\frac{S_w-S_{wi}}{1-S_{wi}-S_{or}})^{m_w}

        k_{rn} = k_{rn, max}(\frac{1-S_w-S_{nr}}{1-S_{wi}-S_{nr}})^{m_n}

    Where

    +-------------------------+-----------------------------------------------+
    |  :math:`m_i`            | Corey Exponents                               |
    +-------------------------+-----------------------------------------------+
    |  :math:`k_r`            | Relative permeability                         |
    +-------------------------+-----------------------------------------------+
    |  :math:`w,n`            | Wetting and NonWetting indices respectively   |
    +-------------------------+-----------------------------------------------+
    |  :math:`S_{wi}`         | Irreducible wetting saturation                |
    +-------------------------+-----------------------------------------------+
    |  :math:`S_{nr}`         | Residual non-wetting saturation               |
    +-------------------------+-----------------------------------------------+
    |  :math:`k_{rn, max}`    | End point relative permeability to oil        |
    |                         | normalized to absolute plug air permeability  |
    +-------------------------+-----------------------------------------------+
    |  :math:`k_{rw, max}`    | End point relative permeability to water      |
    |                         | normalized to absolute plug air permeability  |
    +-------------------------+-----------------------------------------------+

    \endrst

SourceFiles
    krBrooksCorey.C

\*---------------------------------------------------------------------------*/

#ifndef krBrooksCorey_H
#define krBrooksCorey_H

#include "phase.H"
#include "relativePermeabilityModelBase.H"
#include "relativePermeabilityModelBasesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------*\
        Class krBrooksCorey Declaration
\*---------------------------------------------------------------------------*/

class krBrooksCorey
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

        // The other phase
        word theOtherPhase_;

        //- Const-ref to the canonical saturation
        const volScalarField& S_;

        //- Irreducible wetting saturation
        volScalarField Swi_;

        //- Residual nonwetting saturation
        volScalarField Snr_;

        //- Model exponent to wetting & non-wetting phases
        volScalarField mw_, mn_;

        //- End-point permeabilities
        volScalarField krwMax_,krnMax_;

public:

    //- Runtime type information
    TypeName("BrooksCorey");

    // Constructors

        //- Construct from components
        krBrooksCorey
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

        //- Destructor
        ~krBrooksCorey() {}

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

        //- Correct the relative permeabilities
        virtual void correct();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
