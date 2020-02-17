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

#include "zeroGradientFvPatchField.H"

namespace Foam 
{

defineTypeNameAndDebug(wellModel, 0);
defineRunTimeSelectionTable(wellModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class KType, int nPhases>
wellBase<KType, nPhases>::wellBase
(
    const word& name,
    const dictionary& wellProperties,
    const fvMesh& mesh,
    const wellModel& corrector
)
    :
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            corrector,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    driveHandling(nPhases+int(cmpFactorial(nPhases)/2)+1),
    name_(name),
    wellProperties_(wellProperties),
    mesh_(mesh),
    corrector_(corrector),
    wellDict_(wellProperties.subDict(name)),
    operation_
    (
        wordToOperationHandling(wellDict_.lookup("operationMode"))
    ),
    iPhase_
    (
        (operation_ == INJE)
        ? wellDict_.lookupOrDefault<word>("injectionPhase","water")
        : "noPhase"
    ),
    radius_
    (
        dimensionedScalar(name_+".radius", dimLength, readScalar(wellDict_.lookup("radius")))
    ),
    perfos_(),
    wellSet_(mesh, name+"Set", 0),
    phases_
    (
        wellDict_.found("phases") 
        ? wellDict_.lookup("phases")
        : mesh.objectRegistry::names<phase>()
    ),
    source_(),
    tV_(0.0),
    bhp_(name_+".bhp", dimPressure, 0.0),
    tRate_(name_+".targetRate", dimVolume/dimTime, 0.0),
    driveSeries_(),
    timeForDt_(VGREAT)
{
    // Read Well perforations
    // Updates perfos_ and wellSet_
    readPerforations();

    // Update total cells volume
    cellsVolume();

    // Update imposed drives
    readImposedDrives();

    // Initiate Drive HashTable
    driveHandling.insert("BHP", 0);
    forAll(phases_, phi)
    {
        driveHandling.insert(phases_[phi]+".rate", 0);
    }
    List<wordList> phaseCombs = findBinaryCombinations<word>(phases_);
    forAll(phaseCombs, pc)
    {
        driveHandling.insert
        (
            phaseCombs[pc][0]+word(".")+phaseCombs[pc][1]+word(".ratio"), 0
        );
        driveHandling.insert
        (
            phaseCombs[pc][1]+word(".")+phaseCombs[pc][0]+word(".ratio"), 0
        );
    }
    driveHandling.insert("totalRate", 0);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class KType, int nPhases>
autoPtr<wellBase<KType, nPhases> >
wellBase<KType, nPhases>::New
(
    const word& name,
    const dictionary& wellProperties,
    const fvMesh& mesh,
    const wellModel& corrector
)
{
    // Get the name from the dictionary.
    const word modelType(wellProperties.subDict(name).lookup("wellType"));

    // Get the RTS Table via the global object.
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);
    // If the constructor pointer is not found in the table.
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "wellBase::New(const dictionary&)"
        )   << "Unknown well type "
            << modelType << nl << nl
            << "Valid wells are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr< wellBase<KType, nPhases> > 
        (cstrIter()(name, wellProperties, mesh, corrector));
}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

template<class KType, int nPhases>
word wellBase<KType, nPhases>::operationHandlingToWord
(
    const operationHandling& op
) const
{
    word enumName("production");
    switch (op)
    {
        case wellBase<KType, nPhases>::INJE:
        {
            enumName = "injection";
            break;
        }
        case wellBase<KType, nPhases>::PROD:
        {
            enumName = "production";
            break;
        }
    }
    return enumName;
}

template<class KType, int nPhases>
void Foam::wellBase<KType, nPhases>::correct()
{
    preCorrect();
    corrector_(name_);
    postCorrect();
}

template<class KType, int nPhases>
typename wellBase<KType, nPhases>::operationHandling 
Foam::wellBase<KType, nPhases>::wordToOperationHandling
(
    const word& op
) const
{
    if (op == "production")
    {
        return wellBase<KType, nPhases>::PROD;
    }
    else if (op == "injection") 
    {
        return wellBase<KType, nPhases>::INJE;
    }
    else 
    {
        WarningIn
        (
            "Foam::wellBase<KType, nPhases>::wordToOperationHandling( const word& op)"
        )   << "Bad well operation mode specifier " << op << ", using 'production'" << endl;
        return wellBase<KType, nPhases>::PROD;
    }
}


template<class KType, int nPhases>
void wellBase<KType, nPhases>::readPerforations()
{
    Info << "Constructing cells for well: " << name_ << nl;

    // Temp-Store perforation entries in dict
    const PtrList<entry> perfsInfo 
    (
        wellDict_.lookup("perforations")
    );

    // Reshape perforations list
    perfos_.setSize(perfsInfo.size());

    // Construct cells from given topoSetSources
    forAll(perfos_, perfi)
    {
        // Select a perforation interval
        const entry& perfInfo = perfsInfo[perfi];

        // Require that the perforation interval is a valid dict
        if(!perfInfo.isDict())
        {
            FatalIOErrorIn(__PRETTY_FUNCTION__, wellDict_)
                << "Entry " << perfInfo << " in wells section is not a"
                << " valid dictionary." << exit(FatalIOError);
        }

        // Set the pointer to the requested topoSetSource
        perfos_.set(
            perfi,
            topoSetSource::New(perfInfo.keyword(), mesh_, perfInfo.dict())
        );

        // Include selected cells in the well's cell set
        perfos_[perfi].applyToSet(topoSetSource::ADD, wellSet_);
    }

    // Write well set to disk
    if (debug) wellSet_.write();
}


template<class KType, int nPhases>
void wellBase<KType, nPhases>::readImposedDrives()
{
    Info << "\tReading imposed drives for well: " << name_ << nl;

    // Temp-Store drives entries in dict
    const PtrList<entry> drivesInfo 
    (
        wellDict_.lookup("imposedDrives")
    );

    // Reshape drive series list
    driveSeries_.setSize(nPhases-1);

    // Fail if supplied drives are unacceptable
    if (driveSeries_.size() != drivesInfo.size())
        FatalIOErrorIn(__PRETTY_FUNCTION__, wellDict_)
            << "Expected " << driveSeries_.size()
            << " imposedDrives But " << drivesInfo.size()
            << " encountered."
            << exit(FatalError);

    // Construct imposed drives list
    forAll(driveSeries_, di)
    {
        // Select a drive
        const entry& driveInfo = drivesInfo[di];

        // Require that the drive is a valid dict
        if(!driveInfo.isDict())
        {
            FatalIOErrorIn(__PRETTY_FUNCTION__, wellDict_)
                << "Entry " << driveInfo << " in wells section is not a"
                << " valid dictionary." << exit(FatalIOError);
        }

        // Set the requested drive
        driveSeries_[di] = interpolationTable<scalar>(driveInfo.dict());

        //- Activate drive in the hashTable
        driveHandling.set(driveInfo.keyword(), di+1);
    }

}









Foam::wellModel::wellModel
(
    const word& name,
    const dictionary& wellProperties,
    const fvMesh& mesh
)
:
    objectRegistry
    (
        IOobject
        (
            name,
            mesh.time()
        )
    ),
    name_(name),
    wellProperties_(wellProperties),
    mesh_(mesh),
    p_(mesh.lookupObject<volScalarField>("p")),
    source_(p_, dimless/dimTime)
{
}

Foam::autoPtr<Foam::wellModel> Foam::wellModel::New
(
    const word& name,
    const dictionary& wellProperties,
    const fvMesh& mesh
)
{
    const word wellManagement(wellProperties.lookup("wellModel"));

    Info<< "Selecting wellModel : " << wellManagement << "\n" << endl;

    dictionaryConstructorTable::iterator cstrIter =
      dictionaryConstructorTablePtr_->find(wellManagement);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
      {
          FatalErrorIn("Foam::autoPtr<Foam::wellModel> Foam::wellModel::New\
                        (\
                            const word& name,\
                            const dictionary& wellProperties,\
                            const fvMesh& mesh,\
                            const volScalarfield& p\
                        )" )
              << "Unknown wellManagement model "
              << wellManagement << nl << nl
              << "Valid wellManaegement models are : " << endl
              << dictionaryConstructorTablePtr_->sortedToc()
              << exit(FatalError);
      }

    return autoPtr<wellModel>
      (cstrIter()(name, wellProperties, mesh));
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
