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

Namespace
    Foam::relativePermeabilityModels

Description
    A namespace for relative permeability models.

Class
    Foam::relativePermeabilityModelBase

Description
    An abstract class template for relative permeability models.

SourceFiles
    relativePermeabilityModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef relativePermeabilityModelBase_H
#define relativePermeabilityModelBase_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "phase.H"
#include "volFields.H"
#include "volFieldsFwd.H"
#include "krPcCompileCalculations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class relativePermeabilityModelBase Declaration
\*---------------------------------------------------------------------------*/

template<int nPhases>
class relativePermeabilityModelBase
: public regIOobject
{

protected:

    // Protected data

        //- Model Name
        word name_;

        //- Const-ref to Hosting dict
        const dictionary& transportProperties_;

        //- Const-ref to mesh
        const fvMesh& mesh_;

        //- Phase Relative permeabilities
        volVectorNField<numberOfKrFields(nPhases)> kr_;

        //- Phase Relative permeability derivatives
        volVectorNField<numberOfKrFields(nPhases)> dkrdS_;

        //- Disallow copy construct
        relativePermeabilityModelBase(const relativePermeabilityModelBase&);

        //- Disallow default bitwise assignment
        void operator=(const relativePermeabilityModelBase&);

public:

    // Static data members
    TypeName("relativePermeabilityModelBase");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable ( 	  	
        autoPtr,
		relativePermeabilityModelBase,
		dictionary,
		(
        const word& name,
        const dictionary& transportProperties,
        const fvMesh& mesh
        ),
        (name, transportProperties, mesh)
	); 		

    // Constructors

        //- Construct from components
        relativePermeabilityModelBase
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

    // Selectors

        //- Return a reference to the selected Kr model
        static autoPtr<relativePermeabilityModelBase<nPhases> > New
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~relativePermeabilityModelBase() {};

    // Public Member Functions

        //- Return model hosting dictionary
        const dictionary& transportProperties() const
        {
            return transportProperties_;
        }

        //- Return relative permeability component
        tmp<volScalarField> kr(int phaseIndex) const
        {
            return kr_.component(phaseIndex);
        }

        //- Return relative permeability of cell
        scalar kr(int phaseIndex, int cellID) const
        {
            return kr_[phaseIndex][cellID];
        }

        //- Return relative permeability derivative
        tmp<volScalarField> dkrdS(int phaseIndex) const
        {
            return dkrdS_.component(phaseIndex);
        }

        //- Return relative permeability derivative of cell
        scalar dkrdS(int phaseIndex, int cellID) const
        {
            return dkrdS_[phaseIndex][cellID];
        }

        //- Return Kr table as a hole (For debugging)
        const volVectorNField<numberOfKrFields(nPhases)>& krTable() const
        {
            return kr_;
        }

        //- Return Kr table as a hole (For debugging)
        const volVectorNField<numberOfKrFields(nPhases)>& dkrdSTable() const
        {
            return dkrdS_;
        }

        //- Correct relative permeability fields
        virtual void correct() = 0;

        //- Return associated phase names as a list
        virtual const wordList phaseNames() const = 0;

        //- Return canonical phase name
        virtual const word canonicalPhase() const = 0;

        //- Return phase index
        virtual label phaseIndex(const word& phaseName) const = 0;

        //- Write Data to Ostream
        virtual bool writeData(Ostream&) const
        {
            notImplemented(__PRETTY_FUNCTION__);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Macro for RTS-enabled class initiated from this class template
#define makeRelativePermeabilityModelBase(nPhases)                             \
                                                                               \
defineNamedTemplateTypeNameAndDebug(base##nPhases##PhasesKrModel, 0);          \
defineTemplateRunTimeSelectionTable(base##nPhases##PhasesKrModel, dictionary); \

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "relativePermeabilityModelBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
