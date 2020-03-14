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
        base2PhasesKrModel,
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
    const fvMesh& mesh
)
    :
    base2PhasesKrModel(name, transportProperties, mesh),
    coeffsDict_(transportProperties.subDict(typeName+"Coeffs")),
    phaseNames_(
        coeffsDict_.found("phases")
        ? coeffsDict_.lookup("phases")
        : mesh_.objectRegistry::names<phase>()
    ),
    canonicalPhase_(coeffsDict_.lookup("canonicalPhase")),
    theOtherPhase_(phaseNames_[0] != canonicalPhase_ ? phaseNames_[0] : phaseNames_[1]),
    S_(mesh.objectRegistry::lookupObject<volScalarField>(canonicalPhase_+".alpha")),
    Swi_
    (
        IOobject
        (
            canonicalPhase_+".alphaIrr",
            S_.time().timeName()+"/"+typeName,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "Swr", dimless,
            readScalar(coeffsDict_.lookup(canonicalPhase_+".alphaIrr"))
        )
    ),
    Snr_
    (
        IOobject
        (
            theOtherPhase_+".alphaRes",
            S_.time().timeName()+"/"+typeName,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "Snr", dimless,
            readScalar(coeffsDict_.lookup(theOtherPhase_+".alphaRes"))
        )
    ),
    mw_
    (
        IOobject
        (
            canonicalPhase_+".m",
            S_.time().timeName()+"/"+typeName,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            canonicalPhase_+".m", dimless,
            readScalar(coeffsDict_.lookup(canonicalPhase_+".m"))
        )
    ),
    mn_
    (
        IOobject
        (
            theOtherPhase_+".m",
            S_.time().timeName()+"/"+typeName,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            theOtherPhase_+".m", dimless,
            readScalar(coeffsDict_.lookup(theOtherPhase_+".m"))
        )
    ),
    krwMax_
    (
        IOobject
        (
            canonicalPhase_+".krMax",
            S_.time().timeName()+"/"+typeName,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            canonicalPhase_+".krMax", dimless,
            readScalar(coeffsDict_.lookup(canonicalPhase_+".krMax"))
        )
    ),
    krnMax_
    (
        IOobject
        (
            theOtherPhase_+".krMax",
            S_.time().timeName()+"/"+typeName,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            theOtherPhase_+".krMax", dimless,
            readScalar(coeffsDict_.lookup(theOtherPhase_+".krMax"))
        )
    )
{
    // Input Checks
    if (gMin(mn_) <= 0)
    {
        FatalErrorIn(__PRETTY_FUNCTION__)
            << "Corey coefficient m is equal or less than 0" 
            << exit(FatalError);
    }
    if (gMin(mw_) <= 0)
    {
        FatalErrorIn(__PRETTY_FUNCTION__)
            << "Corey coefficient m is equal or less than 0" 
            << exit(FatalError);
    }

    if (gMin(Swi_) < 0 || gMax(Swi_) > 1 )
    {
        FatalErrorIn(__PRETTY_FUNCTION__)
            << "Unexpected saturation " << Swi_.name() << "value" 
            << exit(FatalError);
    }

    if (gMin(Snr_) < 0 || gMax(Snr_) > 1 )
    {
        FatalErrorIn(__PRETTY_FUNCTION__)
            << "Unexpected saturation " << Snr_.name() << "value" 
            << exit(FatalError);
    }

    if (phaseNames_.size() != 2)
    {
        FatalErrorIn(__PRETTY_FUNCTION__)
            << typeName << " Model got unexpected number of phases."
            << " Please set phaseNames keyword in " 
            << coeffsDict_.name() << " dictionary" << endl;
    }
}

// * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * * //

void Foam::relativePermeabilityModels::krBrooksCorey::correct()
{
    forAll(mesh_.C(), celli)
    {
        if( S_[celli] < Swi_[celli] )
        {
            kr_[celli][0] = VSMALL;
            kr_[celli][1] = krnMax_[celli];
            dkrdS_[celli][0] = 0;
            dkrdS_[celli][1] = 0;
        } else if( S_[celli] > (1-Snr_[celli]) ) {
            kr_[celli][0] = krwMax_[celli];
            kr_[celli][1] = VSMALL;
            dkrdS_[celli][0] = 0;
            dkrdS_[celli][1] = 0;
        } else {
            // Calculate Krs
            dimensionedScalar SneUpper = 1-S_.internalField()[celli]-Snr_.internalField()[celli];
            dimensionedScalar SweLower = 1-Swi_.internalField()[celli]-Snr_.internalField()[celli];
            dimensionedScalar SweUpper = S_.internalField()[celli]-Swi_.internalField()[celli];
            kr_.internalField()[celli][0] = 
                krwMax_.internalField()[celli] * pow(SweUpper/SweLower, mw_.internalField()[celli]).value();
            kr_.internalField()[celli][1] = 
                krnMax_.internalField()[celli] * pow(SneUpper/SweLower, mn_.internalField()[celli]).value();

            // Calculate Kr Derivatives irt S_
            if (SneUpper.value() != 0 && SweUpper.value() != 0)
            {
                dkrdS_.internalField()[celli][0] = mw_[celli]*kr_.internalField()[celli][0]/SweUpper.value();
                dkrdS_.internalField()[celli][1] = -mn_[celli]*kr_.internalField()[celli][1]/SneUpper.value();
            } else {
                dkrdS_.internalField()[celli][0] =
                    mw_[celli]*krwMax_.internalField()[celli] * pow(SweUpper/SweLower, mw_[celli]-1).value()/SweLower.value();
                dkrdS_.internalField()[celli][1] =
                    -mn_[celli]*krnMax_.internalField()[celli] * pow(SneUpper/SweLower, mn_[celli]-1).value()/SweLower.value();
            }
        }
    }

    // It's unbelievable how important these lines are !!!!
    kr_.correctBoundaryConditions();
    dkrdS_.correctBoundaryConditions();

}

// ************************************************************************* //
