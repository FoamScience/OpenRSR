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

#include "blackOilFlash.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace bubblePointModels
{
defineTypeNameAndDebug(blackOilFlash, 0);

addToRunTimeSelectionTable
(
    bubblePointModel,
    blackOilFlash,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubblePointModels::blackOilFlash::blackOilFlash
(
    const word& name,
    const dictionary& transportProperties,
    const volScalarField& p
)
    :
    bubblePointModel(name, transportProperties, p),
    PbCoeffs_(transportProperties.subDict(typeName+"Coeffs")),
    Ko_
    (
        IOobject
        (
            "Ko",
            p_.time().timeName(),
            p_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        p_.mesh(),
        dimensionedScalar("", dimless, 0.0)
    ),
    Kg_
    (
        IOobject
        (
            "Kg",
            p_.time().timeName(),
            p_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        p_.mesh(),
        dimensionedScalar("", dimless, 0.0)
    ),
    RsSeries_(PbCoeffs_),
    rhoMost_("rhoMost", dimDensity/dimMoles, readScalar(PbCoeffs_.lookup("rhoMo"))),
    rhoMgst_("rhoMgst", dimDensity/dimMoles, readScalar(PbCoeffs_.lookup("rhoMg"))),
    Ro_
    (
        IOobject
        (
            "Ro",
            p_.time().timeName(),
            p_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        p_.mesh(),
        dimensionedScalar("", dimless, 0.0)
    ),
    Rg_
    (
        IOobject
        (
            "Rg",
            p_.time().timeName(),
            p_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        p_.mesh(),
        dimensionedScalar("", dimless, 0.0)
    ),
    ro_
    (
        IOobject
        (
            "ro",
            p_.time().timeName(),
            p_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        p_.mesh(),
        dimensionedScalar("", dimless, 0.0)
    ),
    rg_
    (
        IOobject
        (
            "rg",
            p_.time().timeName(),
            p_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        p_.mesh(),
        dimensionedScalar("", dimless, 0.0)
    ),
    xo_("xo", ro_),
    yo_("yo", ro_),
    xg_("xg", ro_),
    yg_("yg", ro_)
{
}

// * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * * //

void Foam::bubblePointModels::blackOilFlash::correct()
{
    // First read Rg & Ro from table
    forAll(p_.internalField(), celli)
    {
        Ro_[celli] = RsSeries_(p_.internalField()[celli])[0];
        Rg_[celli] = RsSeries_(p_.internalField()[celli])[1];
    }

    // Calculate molar Rs
    ro_ == Ro_*(rhoMost_/rhoMgst_);
    rg_ == Rg_*(rhoMgst_/rhoMost_);

    // Calculate phase compositions
    // DISCUSSION: Is this absolutely necessary?? Yes it is.
    xo_ = 1.0/(1.0+rg_);
    xg_ = rg_*xo_;
    yg_ = 1.0/(1.0+ro_);
    yo_ = ro_*yg_;

    // Caculate K-values
    Ko_ = (1.0+rg_)/(1.0+1.0/ro_);
    Kg_ = (1.0+1.0/rg_)/(1.0+ro_);

    // Ko*zo+Kg*zg< 1 ---> undersaturated
    // Ko*zo+Kg*zg>=1 ---> saturated
    state_ = sign(Ko_*(xo_+yo_)-Kg_*(xg_+yg_)-1);

    // Fail if free gas is generated
    // DISCUSSION: Why fail?? should I just use time.writeAndEnd()
    //             to conculde the simulation normally??
    forAll(state_.internalField(), celli)
    {
        if (state_[celli] >= 0)
        {
            if (debug) p_.write();
            FatalErrorIn("void Foam::bubblePointModels::constantPb::correct()")
                << "Free gas present at cell" << celli << "." << nl
                << "Check the SaturationState field ..."
                << exit(FatalError);
        }
    }
}

// ************************************************************************* //
