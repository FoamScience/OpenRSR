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

\*---------------------------------------------------------------------------*/

#include "constantPb.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace bubblePointModels
{
defineTypeNameAndDebug(constantPb, 0);

addToRunTimeSelectionTable
(
    bubblePointModel,
    constantPb,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubblePointModels::constantPb::constantPb
(
    const word& name,
    const dictionary& transportProperties,
    const volScalarField& p
)
    :
    bubblePointModel(name, transportProperties, p),
    PbCoeffs_(transportProperties.subDict(typeName+"Coeffs")),
    Pb_("Pb", dimPressure, readScalar(PbCoeffs_.lookup("Pb")))
{
}

// * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * * //

void Foam::bubblePointModels::constantPb::correct()
{
    // State is  1 if Pb>=p  (saturated)
    //       or -1 if Pb< p  (undersaturated)
    state_ = sign(Pb_-p_); 

    forAll(state_.internalField(), celli)
    {
        if (state_[celli] > 0)
        {
            if (debug) p_.write();
            state_.write();
            Info << state_ << endl;
            FatalErrorIn("void Foam::bubblePointModels::constantPb::correct()")
                << "Free gas present at cell " << celli << "." << nl
                << "Check the SaturationState field in Time = "
                << p_.mesh().time().timeName()
                << exit(FatalError);
        }
    }
}

// ************************************************************************* //
