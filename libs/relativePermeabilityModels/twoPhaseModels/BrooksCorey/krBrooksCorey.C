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

#include "krBrooksCorey.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{
defineTypeNameAndDebug(krBrooksCorey, 0);

addToRunTimeSelectionTable
(
    relativePermeabilityModel,
    krBrooksCorey,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModels::krBrooksCorey::krBrooksCorey
(
    const word& name,
    const dictionary& transportProperties,
    const phase& phase1,
    const phase& phase2
)
    :
    relativePermeabilityModel(name, transportProperties, phase1, phase2),
    coeffsDict_(transportProperties.subDict(typeName+"Coeffs")),
    Swc_
    (
        IOobject
        (
            phase1_.name()+".alphaMin",
            S_.time().timeName()+"/"+typeName,
            S_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S_.mesh(),
        dimensionedScalar
        (
            "Swc", dimless,
            readScalar(coeffsDict_.lookup(phase1.name()+".alphaMin"))
        )
    ),
    Sor_
    (
        IOobject
        (
            phase2_.name()+".alphaMin",
            S_.time().timeName()+"/"+typeName,
            S_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S_.mesh(),
        dimensionedScalar
        (
            "Sor", dimless,
            readScalar(coeffsDict_.lookup(phase2.name()+".alphaMin"))
        )
    ),
    n_
    (
        IOobject
        (
            "n",
            S_.time().timeName()+"/"+typeName,
            S_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S_.mesh(),
        dimensionedScalar
        (
            "n", dimless,
            readScalar(coeffsDict_.lookup("n"))
        )
    ),
    Se_((S_-Swc_)/(1.0-Sor_-Swc_)),
    kr1Max_
    (
        IOobject
        (
            phase1_.name()+".krMax",
            S_.time().timeName()+"/"+typeName,
            S_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S_.mesh(),
        dimensionedScalar
        (
            "kr1Max", dimless,
            readScalar(coeffsDict_.lookup(phase1.name()+".krMax"))
        )
    ),
    kr2Max_
    (
        IOobject
        (
            phase2_.name()+".krMax",
            S_.time().timeName()+"/"+typeName,
            S_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S_.mesh(),
        dimensionedScalar
        (
            "kr2Max", dimless,
            readScalar(coeffsDict_.lookup(phase2.name()+".krMax"))
        )
    )
{

    // Input ckecks
    if (gMin(n_) <= 0)
    {
        FatalErrorIn(
            "Foam::relativePermeabilityModels::krBrooksCorey::krBrooksCorey\
            (\
                const word& name,\
                const dictionary& transportProperties,\
                const phase& phase1,\
                const phase& phase2\
            ) ")
            << "Relative permeability coefficient n equal or less than 0" 
            << exit(FatalError);
    }

    if (gMin(Swc_) < 0 || gMax(Swc_) > 1 )
    {
        FatalErrorIn(
            "Foam::relativePermeabilityModels::krBrooksCorey::krBrooksCorey\
            (\
                const word& name,\
                const dictionary& transportProperties,\
                const phase& phase1,\
                const phase& phase2\
            ) ")
            << "Unexpected saturation " << Swc_.name() << "value" 
            << exit(FatalError);
    }

    if (gMin(Sor_) < 0 || gMax(Sor_) > 1 )
    {
        FatalErrorIn(
            "Foam::relativePermeabilityModels::krBrooksCorey::krBrooksCorey\
            (\
                const word& name,\
                const dictionary& transportProperties,\
                const phase& phase1,\
                const phase& phase2\
            ) ")
            << "Unexpected saturation " << Sor_.name() << "value" 
            << exit(FatalError);
    }
    
    // BUG FIX: saturation going under critical value
    // Bound normalized saturation to [0.0,1.0] inclusive
    //forAll(Se_.internalField(), celli)
    //{
    //    if (Se_[celli] <= 0.0) Se_[celli] = 0.0;
    //    if (Se_[celli] >= 1.0) Se_[celli] = 1.0;
    //}
}

// * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * * //

void Foam::relativePermeabilityModels::krBrooksCorey::correct()
{
    if (debug)
    {
        forAll(S_.internalField(), celli)
        {
            if (S_[celli] < Swc_[celli])
                WarningIn
                (
                    "void Foam::relativePermeabilityModels::krBrooksCorey::correct()"
                )
                << "Underflow, " << phase1_.name() 
                << " saturation falls bellow minimal saturation at cell "
                << celli << endl;
            if (S_[celli] > 1.0-Sor_[celli])
                WarningIn
                (
                    "void Foam::relativePermeabilityModels::krBrooksCorey::correct()"
                )
                << "Overflow, " << phase1_.name() 
                << " saturation exceeds maximal saturation at cell "
                << celli << endl;
        }
    }

    // Calculate with BC overwriting
    // Otherwise use correctBoundaryConditions()
    Se_  == (S_-Swc_)/(1.0-Sor_-Swc_);

    kr1_ == kr1Max_ * pow(Se_,n_);
    kr2_ == kr2Max_ * pow((1.0-Se_),n_);
    dkr1dS_ ==  kr1Max_*n_*pow(Se_,(n_-1))/(1.0-Sor_-Swc_);	
    dkr2dS_ == -kr2Max_*n_*pow((1.0-Se_),(n_-1))/(1.0-Sor_-Swc_);

    // Fix the internalField
    forAll(Se_.internalField(), celli)
    {
        if(Se_[celli] < 0.0 or Se_[celli] > 1.0)
        {
            dkr1dS_[celli] = 0.0;
            dkr2dS_[celli] = 0.0;
            kr1_[celli] = 0.0;
            kr2_[celli] = kr2Max_[celli];
            if(Se_[celli] > 1.0)
            {
                kr1_[celli] = kr1Max_[celli];
                kr2_[celli] = 0.0;
            }
        }
    }

    // Fix the boundaryField
    forAll(Se_.boundaryField(), patchi)
    {
        forAll(Se_.boundaryField()[patchi], facei)
        {
            scalar bSe = Se_.boundaryField()[patchi][facei];
            if (bSe < 0.0 or bSe > 1.0)
            {
                dkr1dS_.boundaryField()[patchi][facei] = 0.0;
                dkr2dS_.boundaryField()[patchi][facei] = 0.0;
                kr1_.boundaryField()[patchi][facei] = 0.0;
                kr2_.boundaryField()[patchi][facei] = kr2Max_.boundaryField()[patchi][facei];
                if(bSe > 1.0)
                {
                    kr1_.boundaryField()[patchi][facei] = kr1Max_.boundaryField()[patchi][facei];
                    kr2_.boundaryField()[patchi][facei] = 0.0;
                }

            }
        }
    }

    // It's unbelievable how important these lines are !!!!
    kr1_.correctBoundaryConditions();
    kr2_.correctBoundaryConditions();
    dkr1dS_.correctBoundaryConditions();
    dkr2dS_.correctBoundaryConditions();


}

// ************************************************************************* //
