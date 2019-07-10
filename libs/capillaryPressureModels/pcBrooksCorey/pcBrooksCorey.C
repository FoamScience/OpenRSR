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

#include "pcBrooksCorey.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  namespace capillaryPressureModels
  {
    defineTypeNameAndDebug(pcBrooksCorey, 0);

    addToRunTimeSelectionTable
    (
     capillaryPressureModel,
     pcBrooksCorey,
     dictionary
     );
  }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillaryPressureModels::pcBrooksCorey::pcBrooksCorey
(
    const word& name,
    const dictionary& transportProperties,
    const phase& wettingPhase,
    const phase& nonWettingPhase
)
:
capillaryPressureModel(name, transportProperties, wettingPhase, nonWettingPhase),	
coeffsDict_(transportProperties.subDict(typeName+"Coeffs")),
pcSmin_
(
    IOobject
    (
        Sw_.name()+"PcMin",
        Sw_.time().timeName()/typeName,
        Sw_.db(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    Sw_.mesh(),
    dimensionedScalar(
        "", 
        dimless, 
        coeffsDict_.lookupOrDefault<scalar>(wettingPhase.name()+".pcSmin",0.0)
    )
),
pcSmax_
(
    IOobject
    (
        Sw_.name()+"PcMax",
        Sw_.time().timeName()/typeName,
        Sw_.db(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    Sw_.mesh(),
    dimensionedScalar(
        "", 
        dimless, 
        coeffsDict_.lookupOrDefault<scalar>(wettingPhase.name()+".pcSmax",1.0)
    )
),
pc0_
(
    IOobject
    (
        "pc0",
        Sw_.time().timeName()/typeName,
        Sw_.db(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    Sw_.mesh(),
    dimensionedScalar(
        "", 
        dimPressure,
        coeffsDict_.lookupOrDefault<scalar>("pc0",0.0)
    )
),
alpha_
(
    IOobject
    (
        "alpha",
        Sw_.time().timeName()/typeName,
        Sw_.db(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    Sw_.mesh(),
    dimensionedScalar(
        "", 
        dimless, 
        coeffsDict_.lookupOrDefault<scalar>("pc.alpha",1)
    )
),
Se_((Sw_- pcSmin_)/(pcSmax_-pcSmin_))
{
    // Input checks
    if (gMin(alpha_) <= 0) 
        FatalErrorIn("Foam::capillaryPressureModels::pcBrooksCorey::pcBrooksCorey\
                    (\
                        const word& name,\
                        const dictionary& transportProperties,\
                        const phase& wettingPhase,\
                        const phase& nonWettingPhase\
                    )")
            << "alpha < 0 in BrooksCorey capillary presssure model." 
            << exit(FatalError);

    if (gMin(pcSmin_) < 0 || gMax(pcSmin_) > 1 )
    {
        FatalErrorIn("Foam::capillaryPressureModels::pcBrooksCorey::pcBrooksCorey\
                    (\
                        const word& name,\
                        const dictionary& transportProperties,\
                        const phase& wettingPhase,\
                        const phase& nonWettingPhase\
                    )")
            << "Unexpected saturation " << pcSmin_.name() << " value" 
            << exit(FatalError);
    }

    if (gMin(pcSmax_) <= 0 || gMax(pcSmax_) >= 1 )
    {
        FatalErrorIn("Foam::capillaryPressureModels::pcBrooksCorey::pcBrooksCorey\
                    (\
                        const word& name,\
                        const dictionary& transportProperties,\
                        const phase& wettingPhase,\
                        const phase& nonWettingPhase\
                    )")
            << "Unexpected saturation " << pcSmax_.name() << " value" 
            << exit(FatalError);
    }

}

void Foam::capillaryPressureModels::pcBrooksCorey::correct()
{
    Se_ == (Sw_-pcSmin_)/(pcSmax_ - pcSmin_);
    pc_ = pc0_ * pow(Se_,alpha_);
    dpcdS_= -alpha_*pc0_*(pow(Se_,alpha_-1))/(pcSmax_-pcSmin_);

    //// Fix the internalField
    //forAll(Se_.internalField(), celli)
    //{
    //    if(Se_[celli] < 0.0)
    //    {
    //        dpcdS_[celli] = 0.0;
    //        pc_[celli] = 0.0;
    //    }
    //    if(Se_[celli] > 1.0)
    //    {
    //        dpcdS_[celli] = 0.0;
    //        pc_[celli] = pc0_[celli];
    //    }
    //}
    //// Fix the boundaryField
    //// DISCUSSION: Does this work through Processor boundaries
    //forAll(Se_.boundaryField(), patchi)
    //{
    //    forAll(Se_.boundaryField()[patchi], facei)
    //    {
    //        scalar bSe = Se_.boundaryField()[patchi][facei];
    //        if(bSe < 0.0)
    //        {
    //            dpcdS_.boundaryField()[patchi][facei] = 0.0;
    //            pc_.boundaryField()[patchi][facei] = 0.0;
    //        }
    //        if(bSe > 1.0)
    //        {
    //            dpcdS_.boundaryField()[patchi][facei] = 0.0;
    //            pc_.boundaryField()[patchi][facei] = pc0_.boundaryField()[patchi][facei];
    //        }
    //    }
    //}
    pc_.correctBoundaryConditions();
    dpcdS_.correctBoundaryConditions();
}
// ************************************************************************* //
