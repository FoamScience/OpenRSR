/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class KType, class PhaseType, int nPhases>
void Foam::wellModelBase<KType, PhaseType, nPhases>::readWells()
{
	notImplemented("void wellModelBase::readWells()");
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class KType, class PhaseType, int nPhases>
Foam::wellModelBase<KType, PhaseType, nPhases>::wellModelBase
(
    const word& name,
    const dictionary& wellModelBaseDict,
    const fvMesh& mesh
)
:
    name_(name),
    wellModelBaseDict_(wellModelBaseDict),
    mesh_(mesh),
    p_(mesh.lookupObject<volScalarField>("p")),
    K_(mesh.lookupObject<KType>("K")),
	phases_(nPhases),
    wells_(),
    sourceMatrices_(nPhases),
    g_
    (
        mesh.objectRegistry::lookupObject<uniformDimensionedVectorField>("g")
    )
{
	// Get phase names 
	// ( Assuming all phases are of the same concrete type 
	// and that type inherits from the "phase" class ).
	const wordList& phaseNames = mesh.objectRegistry::names<phase>();

	// TODO: should check if more than nPhases are requested in wells.
	// ??Potential existence of auxilary phases that don't need well sources??
	
	forAll(phaseNames, pn)
	{
		// Find each phase obj using its name
		phases_[pn].insert
		(
			phaseNames[pn],
		    autoPtr<PhaseType> (mesh.lookupObject<PhaseType>(phaseNames[pn]))
		);
		// Populate source hash table
        sourceMatrices_.insert
		(
			phaseNames[pn], 
			fvScalarMatrix(p_, dimVolume/dimTime)
		);
	}

    // Create well objects and store their pointers to a list
    readWells();
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class KType, class PhaseType, int nPhases>
Foam::autoPtr<Foam::wellModelBase<KType, PhaseType, nPhases> > 
Foam::wellModelBase<KType, PhaseType, nPhases>::New
(
    const word& name,
    const dictionary& wellModelBaseDict,
    const fvMesh& mesh
)
{
    const word wellModelBaseType = wellModelBaseDict.lookup("wellModel");

    Info<< "Selecting wellModel : " << wellModelBaseType << "\n" << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
      dictionaryConstructorTablePtr_->find(wellModelBaseType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
      {
          FatalErrorIn(__PRETTY_FUNCTION__)
              << "Unknown wellModel type "
              << wellModelBaseType << nl << nl
              << "Valid well models are : " << endl
              << dictionaryConstructorTablePtr_->sortedToc()
              << exit(FatalError);
      }

    return autoPtr<wellModelBase<KType, PhaseType, nPhases> >
        (cstrIter()(name, wellModelBaseDict, mesh));
}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class KType, class PhaseType, int nPhases>
const Foam::tmp<Foam::scalarField> 
Foam::wellModelBase<KType, PhaseType, nPhases>::explicitSource
(
	const Foam::word& phaseName
) const
{
	notImplemented("explicitSource(phaseName)");
	return scalarField::zero;
}

// ************************************************************************* //
