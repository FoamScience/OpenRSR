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

#include "well.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(well, 0);
defineRunTimeSelectionTable(well, dictionary);
PtrList<well> well::wellsList;
wordList well::wellsNames;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::well::well
(
    const word& name,
    const dictionary& wellDict,
    const fvMesh& mesh,
    const HCmixture<blackoilPhase>& mixture
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    wellDict_(wellDict),
    mesh_(mesh),
    mixture_(mixture),
    type_(wellDict.lookupOrDefault<word>("type", "verticalWell")),
    mode_
    (
        wordToRateHandling(wellDict.lookup("controlMode"))
    ),
    operation_
    (
        wordToOperationHandling(wellDict.lookup("operationMode"))
    ),
    iPhase_
    (
        (operation_ == INJE)
        ? wellDict.lookupOrDefault<word>("injectionPhase","water")
        : "noPhase"
    ),
    radius_
    (
        dimensionedScalar(name_+".radius", dimLength, readScalar(wellDict_.lookup("radius")))
    ),
    perfos_(),
    wellSet_(mesh, name_+"Set", 0),
    tV_(0.0),
    flowRateData_(wellDict),
    ewSource_(),
    eoSource_(),
    iwSource_(),
    ioSource_(),
    perfosNotRead_(true),
    timeForDt_(VGREAT),
    allowSwitchingModes_(wellDict.lookupOrDefault("allowSwitchingModes",false)),
    switchedToBHP_(false),
    WJw_(mesh.objectRegistry::lookupObject<volScalarField>("p"), dimless/dimTime)
{
}


Foam::well::flowRateModeVars::flowRateModeVars
(
    const dictionary& wellDict
)
:
    flowRateSeries(wellDict.subDict("flowRateData"))
{
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

Foam::word Foam::well::operationHandlingToWord
(
    const operationHandling& op
) const
{
    word enumName("production");

    switch (op)
    {
        case well::INJE:
        {
            enumName = "injection";
            break;
        }
        case well::PROD:
        {
            enumName = "production";
            break;
        }
    }

    return enumName;
}

Foam::word Foam::well::rateHandlingToWord
(
    const rateHandling& rate
) const
{
    word enumName("liquidRate");

    switch (rate)
    {
        case well::LRATE:
        {
            enumName = "liquidRate";
            break;
        }
        case well::BHP:
        {
            enumName = "BHP";
            break;
        }
    }

    return enumName;
}

Foam::well::operationHandling 
Foam::well::wordToOperationHandling
(
    const word& op
) const
{
    if (op == "production")
    {
        return well::PROD;
    }
    else if (op == "injection") 
    {
        return well::INJE;
    }
    else 
    {
        WarningIn
        ( 
            "Foam::well::wordToOperationHandling ( const word& op)"
        )   << "Bad well operation mode specifier " << op << " using 'production'" << endl;
        return well::PROD;
    }
}

Foam::well::rateHandling 
Foam::well::wordToRateHandling
(
    const word& rate
) const
{
    if (rate == "totalFlowRate")
    {
        return well::LRATE;
    }
    else if (rate == "BHP") 
    {
        return well::BHP;
    }
    else 
    {
        WarningIn
        ( 
            "Foam::well::wordToRateHandling ( const word& rate)"
        )   << "Bad well rate control specifier " << rate << " using 'totalFlowRate'" << endl;
        return well::LRATE;
    }
}

Foam::well::operationHandling Foam::well::handleOperationMode
(
    const operationHandling newOp
)
{
    operationHandling prev = operation_;
    operation_ = newOp;
    
    if (debug)
    {
        Info<< "Switching from " << prev << " to " 
            << newOp << " operation mode for well " << name_ 
            << endl;
    }

    return prev;
}

Foam::well::rateHandling Foam::well::handleRateMode
(
    const rateHandling& newRate
)
{
    rateHandling prev = mode_;
    mode_ = newRate;
    
    if (debug)
    {
        Info<< "Switching from " << prev << " to " 
            << newRate << " rate mode for well " << name_ 
            << endl;
    }

    return prev;
}

void Foam::well::readPerforations()
{
    Info << "Constructing cells for well: " << name_ << endl;

    // Store entries
    const PtrList<entry> perfsInfo(
        wellDict_.lookup("perforations")
    );

    // Reshape perfs list
    perfos_.setSize(perfsInfo.size());

    // Construct cells
    forAll(perfos_, perfi)
    {
        const entry& perfInfo = perfsInfo[perfi];

        if(!perfInfo.isDict())
        {
            FatalIOErrorIn("void Foam::well::readPerforations()", wellDict_)
                << "Entry " << perfInfo << " in wells section is not a"
                << " valid dictionary." << exit(FatalIOError);
        }

        perfos_.set(
            perfi,
            topoSetSource::New(perfInfo.keyword(), mesh_, perfInfo.dict())
        );

        perfos_[perfi].applyToSet(topoSetSource::ADD, wellSet_);
    }

    // Write well set to disk
    wellSet_.write();

    // Calculate totalVolume
    // This has to execute once per processor
    forAll(wellSet_.toc(), celli)
    {
        tV_ += mesh_.V()[wellSet_.toc()[celli]];
    }

    Info << "Total Volume for well " << name_ << ":" << tab << tV_ << endl;

    perfosNotRead_ = false;
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::autoPtr<Foam::well> Foam::well::New
(
    const word& name,
    const dictionary& wellDict,
    const fvMesh& mesh,
    const HCmixture<blackoilPhase>& mixture
)
{
    const word wellOrientation(wellDict.lookup("orientation"));

    //Info<< "Selecting well : " << wellOrientation << "\n" << endl;

    dictionaryConstructorTable::iterator cstrIter =
      dictionaryConstructorTablePtr_->find(wellOrientation);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
      {
          FatalErrorIn("Foam::autoPtr<Foam::well> Foam::well::New\
                        (\
                            const word& name,\
                            const dictionary& wellDict,\
                            const fvMesh& mesh\
                        )" )
              << "Unknown well orientation "
              << wellOrientation << nl << nl
              << "Valid well orientations are : " << endl
              << dictionaryConstructorTablePtr_->sortedToc()
              << exit(FatalError);
      }

    return autoPtr<well>
      (cstrIter()(name, wellDict, mesh, mixture));
}

void Foam::well::readWells
(
    const dictionary& wellsProperties,
    const fvMesh& mesh,
    const HCmixture<blackoilPhase>& mixture
)
{
    Info << "Reading Wells" << endl;

    // Store well entries
    const PtrList<entry> wellsInfo(
        wellsProperties.lookup("wells")
    );

    // Reshape wells list
    wellsList.setSize(wellsInfo.size());
    wellsNames.setSize(wellsInfo.size());

    // Construct well ptrs
    forAll(wellsList, welli)
    {
        const entry& wellInfo = wellsInfo[welli];

        if(!wellInfo.isDict())
        {
            FatalIOErrorIn("void Foam::wellModel::readWells()", wellsProperties)
                << "Entry " << wellInfo << " in wells section is not a"
                << " valid dictionary." << exit(FatalIOError);
        }

        wellsList.set(
            welli, 
            well::New(wellInfo.keyword(), wellInfo.dict(), mesh, mixture)
        );
        wellsNames[welli] = wellInfo.keyword();

        if (wellsList[welli].perfosNotRead())
        {
            wellsList[welli].readPerforations(); 
        }
    }
}

// ************************************************************************* //
