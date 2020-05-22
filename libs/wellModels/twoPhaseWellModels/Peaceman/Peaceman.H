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

Class
    Foam::Peaceman

Description
    Peaceman class

SourceFiles
	Peaceman.C

\*---------------------------------------------------------------------------*/

#ifndef Peaceman_H
#define Peaceman_H

#include "wellModelBase.H"
#include "capillaryPressureModelBase.H"
#include "relativePermeabilityModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

namespace twoPhaseWellModels {

/*---------------------------------------------------------------------------*\
                         Class Peaceman Declaration
\*---------------------------------------------------------------------------*/

// Forward declare classes

template <class KType, class PhaseType>
class Peaceman
    : public wellModelBase<KType, PhaseType, 2> {

protected:

    // Protected data

		//- Kr Model name
		const word& krModelName_;

		//- Const-ref to the 2-Phase Kr Model object
		const relativePermeabilityModelBase<2>& krModel_;

		//- Const-ref to the 2-Phase Pc Model object
		const capillaryPressureModelBase<2>& pcModel_;

    // Protected Member Functions and Operators

		//- Calculate Fractional Flow in Well Cells
		scalarList calculateFracFlow(const wellBase<KType, PhaseType, 2>& well) const;

		//- Calculate the ratio: totalRate/rateThroughCell for a well by phase
		scalar calculateTotalRateOverCellRate
		(
			const label& cellID,
			const wellBase<KType, PhaseType, 2>& well,
			const word& phaseName
		) const;

		//- Estimate Well cell sizes
		List<vector> estimateCellSizes(const wellBase<KType, PhaseType, 2>& well) const;

		//- Estimate Well Equivalent radius
		scalarList estimateEquivRadius(const wellBase<KType, PhaseType, 2>& well) const;

		//- Calculate Well Productivity Index
		scalarList calculateWellPI(const wellBase<KType, PhaseType, 2>& well) const;


public:

    //- Runtime type information
    TypeName("Peaceman");

    // Declare Run-Time selection table
    //declareRunTimeSelectionTable
	//(
    //    autoPtr,
    //    Peaceman,
    //    dictionary,
    //    (
    //        const word& name,
    //        const dictionary& PeacemanDict,
    //        const fvMesh& mesh),
    //    (name, PeacemanDict, mesh)
	//);

    // Constructors

		//- Construct from components
    	Peaceman
		(
    	    const word& name,
    	    const dictionary& PeacemanDict,
    	    const fvMesh& mesh
		);

    // Selectors

		//- Return a reference to the selected Peaceman
    	//static autoPtr<Peaceman> New
		//(
    	//    const word& name,
    	//    const dictionary& PeacemanDict,
    	//    const fvMesh& mesh
		//);

    //- Destructor
    virtual ~Peaceman() {};

    // Public Member Functions

    	//- Correct Model Coeffs
    	virtual void correct();

		//- Calculate Well Source for a single well
		virtual void operator()(const word& wellName) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Macro for registering Peaceman models with base well model class 
#define makePeacemanModel(KType, PhaseType)                                   \
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
	      (derived##KType####PhaseType##Peaceman, 0);                         \
    addToRunTimeSelectionTable                                                \
	(                                                                         \
		base2##KType####PhaseType##WellModel,                                 \
		derived##KType####PhaseType##Peaceman,                                \
		Dictionary                                                            \
	);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseWellModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "Peaceman.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
