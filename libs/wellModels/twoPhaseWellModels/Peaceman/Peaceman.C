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

#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include "wellModelBase/wellModelBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //



// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class KType, class PhaseType>
Foam::scalarList Foam::twoPhaseWellModels::Peaceman<KType, PhaseType>
                 ::calculateFracFlow
(
	const Foam::wellBase<KType, PhaseType, 2>& well
) const
{
	notImplemented(__PRETTY_FUNCTION__);
}


template<class KType, class PhaseType>
Foam::scalar Foam::twoPhaseWellModels::Peaceman<KType, PhaseType>
                 ::calculateTotalRateOverCellRate
(
	const Foam::label& cellID,
	const Foam::wellBase<KType, PhaseType, 2>& well,
	const Foam::word& phaseName
) const;
{
	notImplemented(__PRETTY_FUNCTION__);
}


template<class KType, class PhaseType>
List<vector> Foam::twoPhaseWellModels::Peaceman<KType, PhaseType>
                 ::estimateCellSizes
(
	const Foam::wellBase<KType, PhaseType, 2>& well
) const;
{
	notImplemented(__PRETTY_FUNCTION__);
}

template<class KType, class PhaseType>
scalarList Foam::twoPhaseWellModels::Peaceman<KType, PhaseType>
                 ::estimateEquivRadius
(
	const Foam::wellBase<KType, PhaseType, 2>& well
) const;
{
	notImplemented(__PRETTY_FUNCTION__);
}

template<class KType, class PhaseType>
scalarList Foam::twoPhaseWellModels::Peaceman<KType, PhaseType>
                 ::calculateWellPI
(
	const Foam::wellBase<KType, PhaseType, 2>& well
) const;
{
	notImplemented(__PRETTY_FUNCTION__);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class KType, class PhaseType>
Foam::Peaceman<KType, PhaseType>::Peaceman
(
    const word& name,
    const dictionary& PeacemanDict,
    const fvMesh& mesh
)
:
	Foam::wellModelBase<KType, PhaseType, 2>
	(
		name,
		PeacemanDict,
		mesh
	),
	krModelName_
	(
		PeacemanDict.found("names")
		? wellProperties.subDict("names").lookupOrDefault<word>("krModel","krModel")
		: "krModel"
	),
	krModel_
	(
		mesh.objectRegistry::lookupObject<relativePermeabilityModelBase<2> >
		(krModelName_)
	),
	pcModel_
	(
		mesh.objectRegistry::lookupObject<capillaryPressureModelBase<2> >
		(krModelName_)
	),
{
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

template<class KType, class PhaseType>
void Foam::Peaceman<KType, PhaseType>::()
{
}

// ************************************************************************* //
