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
    Foam::capillaryPressureModels::normalizedLogarithmic

Description
    A permeability-independent standard logarithmic mode for capillary
    pressure
    \f[
        P_c = -P_{c0} \log_{10}(S_{w, norm})
    \f]
  
    The normalized saturation is given by:
    \f[
        S_{w, norm} = \frac{S_w - S_{pc, min}}{S_{pc, max}-S_{pc, min}}
    \f]

    Where \f$S_{pc, min} <= S_{w, irr}\f$
    \vartable
        P_c                   | Capillary pressure
        S_w                   | Wetting phase saturation
        S_{pc, norm}          | Normalized capillary saturation
        S_{pc, min}           | Min. capillary saturation
        S_{pc, max}           | Max. cappilary saturation
        \frac{1}{\alpha}      | Pore-size distribution index        
    \endtable

Note
    The derivate \f$\frac{dP_c}{dS_w}\f$ is calculated theoretically.

SourceFiles
    normalizedLogarithmic.C

\*---------------------------------------------------------------------------*/

#ifndef normalizedLogarithmic_H
#define normalizedLogarithmic_H

#include "capillaryPressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace capillaryPressureModels
{

/*---------------------------------------------------------------------------*\
        Class normalizedLogarithmic Declaration
\*---------------------------------------------------------------------------*/

class normalizedLogarithmic
:
    public capillaryPressureModel
{

    //- Configuration dictionary
    dictionary coeffsDict_;

    //- Model fields
    volScalarField pcSmin_;
    volScalarField pcSmax_;
    volScalarField pc0_;

    // effective saturation of phase b
    volScalarField Se_;
        
public:

    //- Runtime type information
    TypeName("normalizedLogarithmic");

    // Constructors

    //- Construct from components
    normalizedLogarithmic
    (
        const word& name,
        const dictionary& transportProperties,
        const phase& wettingPhase,
        const phase& nonWettingPhase
    );

    //- Destructor
    ~normalizedLogarithmic()
        {}

    // Member Functions

    //- Correct the capillary pressure
    void correct();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillaryPressureModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
