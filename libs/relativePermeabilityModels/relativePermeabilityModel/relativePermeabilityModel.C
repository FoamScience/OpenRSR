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

#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(relativePermeabilityModel, 0);
defineRunTimeSelectionTable(relativePermeabilityModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModel::relativePermeabilityModel
(
    const word& name,
    const dictionary& transportProperties,
    const phase& phase1,
    const phase& phase2,
    const word& keyword
)
    :
    name_(name),
    transportProperties_(transportProperties),
    phase1_(phase1),
    phase2_(phase2),
    S_(phase1.alpha()),
    kr1_
    (
        IOobject
        (
            phase1_.name()+".kr",
            S_.time().timeName(),
            S_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S_.mesh(),
        dimensionedScalar("kr", dimless, 0.0)
    ),
    kr2_
    (
        IOobject
        (
            phase2_.name()+".kr",
            S_.time().timeName(),
            S_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S_.mesh(),
        dimensionedScalar("kr", dimless, 0.0)
    ),
    dkr1dS_
    (
        IOobject
        (
            phase1_.name()+".dkrdS",
            S_.time().timeName(),
            S_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S_.mesh(),
        dimensionedScalar("dkrdS", dimless, 0.0)
    ),
    dkr2dS_
    (
        IOobject
        (
            phase2_.name()+".dkrdS",
            S_.time().timeName(),
            S_.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        S_.mesh(),
        dimensionedScalar("dkrdS", dimless, 0.0)
    ),
    wells_(),
    wellKrSwmin_(0), wellKrSnmin_(1),
    wellKrNw_(2), wellKrNo_(2),
    wellKr1Max_(1), wellKr2Max_(1)
{}

Foam::autoPtr<Foam::relativePermeabilityModel> Foam::relativePermeabilityModel::New
(
    const word& name,
    const dictionary& transportProperties,
    const phase& phase1,
    const phase& phase2,
    const word& keyword
)
{
    const word modelType(transportProperties.lookup(keyword));

    Info<< "Selecting relativePermeability model: " 
        << modelType << "\n" << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn(
            "Foam::autoPtr<Foam::relativePermeabilityModel> Foam::relativePermeabilityModel::New\
            (\
                const word& name,\
                const dictionary& transportProperties,\
                const phase& phase1,\
                const phase& phase2\
            ) ")   
            << "Unknown relativePermeabilityModel type "
            << modelType << nl << nl
            << "Valid relativePermeabilityModels are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<relativePermeabilityModel>
        (cstrIter()(name, transportProperties, phase1, phase2));
}


void Foam::relativePermeabilityModel::modifyKrInWells()
{
    // By the time this function is called, wells must be constructed
    // This function should be executed just before correcting BC
    // in all models' correct method
    wells_ = well::wellsNames;
    forAll(wells_, welli)
    {
        //- Const Ref to well object
        const well& currentWell = S_.mesh().lookupObject<well>(wells_[welli]);

        if
        (
             (currentWell.allowSwitchingModes())
                and 
             (not(transportProperties_.found("wellKrModel<"+currentWell.name()+">")))
        )
        {
            // Should we output this only for debugging??
            WarningIn("void Foam::relativePermeabilityModel::modifyKrInWells()")
                << nl << tab
                << "Well " << currentWell.name() << " was set to allow switching "
                << "operation modes but no special Kr model was specified." 
                << nl << tab 
                << "If a segFault occurs, please set a Brooks-Corey Kr coefficients for this well"
                << " in transportPeroperties\n\tor disable switching operation mode:"
                << nl << tab 
                << "wellKrModel<" << currentWell.name() << "> Swmin Somax nw no krwMax kroMax;"
                << nl << endl;
        }

        //- Make changes only if the well is expected to switch operation mode
        if
        (
            (currentWell.allowSwitchingModes())
            and 
            (transportProperties_.found("wellKrModel<"+currentWell.name()+">"))
        )
        {
        transportProperties_.lookup("wellKrModel<"+currentWell.name()+">")
            >> wellKrSwmin_ >> wellKrSnmin_
            >> wellKrNw_ >> wellKrNo_ 
            >> wellKr1Max_ >> wellKr2Max_;

            // The actual modification of kr
            const labelList& cells = currentWell.cellIDs();
            forAll(cells, celli)
            {
                const label cellID = cells[celli];
                // Skip cell if S > Sc
                //if(S_[cellID] > wellKrSnmin_) continue;
                if(kr1_[cellID] > SMALL) continue;

                // Otherwise, correct permeabilites
                scalar Se = (S_[cellID]-wellKrSwmin_)/(1-wellKrSwmin_-wellKrSnmin_);
                kr1_[cellID] = wellKr1Max_*pow(Se, wellKrNw_);
                kr2_[cellID] = wellKr2Max_*pow(1.0-Se, wellKrNo_);
                dkr1dS_[cellID] = wellKr1Max_*wellKrNw_*pow(Se, wellKrNw_-1)/(1-wellKrSwmin_-wellKrSnmin_);
                dkr2dS_[cellID] = -wellKr2Max_*wellKrNo_*pow(Se, wellKrNo_-1)/(1-wellKrSwmin_-wellKrSnmin_);
            }
        }

    }
}
// ************************************************************************* //
