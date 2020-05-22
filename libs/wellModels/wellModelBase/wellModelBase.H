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

Namespace
    Foam::twoPhaseWellModels

Description
    twoPhaseWellModels namespace

Namespace
    Foam::threePhaseWellModels

Description
    threePhaseWellModels namespace

Class
    Foam::wellModelBase

Description
    wellModelBase class

SourceFiles
	wellModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef wellModelBase_H
#define wellModelBase_H

#include "fvMesh.H"
#include "fvMatricesFwd.H"
#include "runTimeSelectionTables.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

/*---------------------------------------------------------------------------*\
                         Class wellModelBase Declaration
\*---------------------------------------------------------------------------*/

// Forward declare classes
class phase;
template<class KType, class PhaseType, int nPhases> class wellBase;

template <class KType, class PhaseType, int nPhases>
class wellModelBase
    : public objectRegistry {

protected:

    // Protected data

		//- Model Name
    	word name_;

    	//- Const-ref to Hosting dict
    	const dictionary& wellModelBaseDict_;

    	//- Const-ref to the mesh where the wells are defined
    	const fvMesh& mesh_;

        //- Const ref to (global/main) pressure field
        const volScalarField& p_;

        //- Const ref to absolute permeability field
        const KType& K_;

        //- Table of Pointers to associated phases
        HashTable<autoPtr<PhaseType> > phases_;

        //- List of Pointers to wells
        PtrList<wellBase<KType, PhaseType, nPhases> > wells_;

        //- Sources matrices
        HashTable<fvScalarMatrix> sourceMatrices_;

        //- Ref to Gravitational field
        const uniformDimensionedVectorField& g_;

    // Protected Member Functions and Operators

    	//- Disallow copy construct
    	wellModelBase(const wellModelBase&);

    	//- Disallow default bitwise assignment
    	void operator=(const wellModelBase&);

		//- Read Wells from object dictionary
		void readWells();

public:

    //- Runtime type information
    TypeName("wellModelBase");

    // Declare Run-Time selection table
    declareRunTimeSelectionTable
	(
        autoPtr,
        wellModelBase,
        dictionary,
        (
            const word& name,
            const dictionary& wellModelBaseDict,
            const fvMesh& mesh),
        (name, wellModelBaseDict, mesh)
	);

    // Constructors

		//- Construct from components
    	wellModelBase
		(
    	    const word& name,
    	    const dictionary& wellModelBaseDict,
    	    const fvMesh& mesh
		);

    // Selectors

		//- Return a reference to the selected wellModelBase
    	static autoPtr<wellModelBase> New
		(
    	    const word& name,
    	    const dictionary& wellModelBaseDict,
    	    const fvMesh& mesh
		);

    //- Destructor
    virtual ~wellModelBase() {};

    // Public Member Functions

		//- Return model name
    	const word& name() const
    	{
    	    return name_;
    	}

    	//- Return model hosting dictionary
    	const dictionary& dict() const
    	{
    	    return wellModelBaseDict_;
    	}

    	//- Return pointer to selected phase
    	const autoPtr<PhaseType> ptrPhases(const word& phaseName) const
    	{
    	    return phases_[phaseName];
    	}

        //- Return Points to wells list
        const PtrList<wellBase<KType, PhaseType, nPhases> >& wells() const {
            return wells_;
        }

        //- Return a ref to model gravitational field
        const uniformDimensionedVectorField& g() const {
            return g_;
        }

        //- Return the object registry
        const objectRegistry& thisDb() const
        {
            return *this;
        }

        //- Return global wells source matrix by phase in implicit settings
		//  i.e, the matrix: A.p = b
        const fvScalarMatrix& source(const word& phaseName) const {
            return sourceMatrices_[phaseName];
        }

        //- Return global wells source by phase in explicit settings as a field
        const tmp<scalarField> explicitSource(const word& phaseName) const;

    	//- Correct Model Coeffs
    	virtual void correct() = 0;

		//- Calculate Well Source for a single well
		virtual void operator()(const word& wellName) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

 // Macro for RTS-enabled class initiated from this class template
#define makeWellModelBase(KType, PhaseType, nPhases)                          \
                                                                              \
	defineNamedTemplateTypeNameAndDebug                                       \
          (base##nPhases####PhaseType####KType##WellModel,0)                  \
    defineTemplateRunTimeSelectionTable                                       \
          (base##nPhases####KType##WellModel, Dictionary);             

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#include "wellModelBase.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
