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
    Foam::wells

Description
    A namespace for well classes and well models.

Class
    Foam::wellBase

Description
    An abstract class template for well objects.

Class
    Foam::wellModel

Description
    An abstract class for well models.

SourceFiles
    wellBase.C

\*---------------------------------------------------------------------------*/

#ifndef wellBase_H
#define wellBase_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "VectorN.H"
#include "volFields.H"
#include "volFieldsFwd.H"
#include "parameterTypes.H"
#include "cellSet.H"
#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class wellBase Declaration
\*---------------------------------------------------------------------------*/

class wellModel;

template<class KType, int nPhases>
class wellBase
    : public regIOobject
{

public:

    // Public data members

        //- Enueration for handling operationModes
        enum operationHandling
        {
            PROD,             /*!< Production well */
            INJE,             /*!< Injection well */
        };

protected:

    // Protected data

        //- Well Name
        word name_;

        //- Const-ref to global web configs
        const dictionary& wellProperties_;

        //- Const-ref to mesh
        const fvMesh& mesh_;

        //- Const-ref to the model driving this well
        const wellModel& corrector_;

        //- Well dictionary
        dictionary wellDict_;

        //- Well operation mode: injection - production
        operationHandling operation_;

        //- Injection phase name
        word iPhase_;

        //- well radius
        dimensionedScalar radius_;

        //- Well Perforations
        PtrList<topoSetSource> perfos_;

        //- Well topo set
        cellSet wellSet_;

        //- Well source
        List<VectorN<scalar, 2*nPhases> > source_;

        //- Well total volume
        scalar  tV_;

       //- TimeStep to catch imposed wellRate/BHP changes
       scalar timeForDt_;

        //- Disallow copy construct
        wellBase(const wellBase&);

        //- Disallow default bitwise assignment
        void operator=(const wellBase&);

public:

    // Static data members
    TypeName("wellBase");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable ( 	  	
        autoPtr,
		wellBase,
		dictionary,
		(
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh,
            const wellModel& corrector
        ),
        (name, wellProperties, mesh, corrector)
	); 		

    // Constructors

        //- Construct from components
        wellBase
        (
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh,
            const wellModel& corrector
        );

    // Selectors

        //- Return a reference to the selected Kr model
        static autoPtr<wellBase<KType, nPhases> > New
        (
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh,
            const wellModel& corrector
        );

    //- Destructor
    virtual ~wellBase() {};

    // Public Member Functions

        //- Return operationHandling as word
        word operationHandlingToWord(const operationHandling& op) const;

        //- Return operationHandling as enumeration
        operationHandling wordToOperationHandling(const word& op) const;

        //- Handle operation mode from enum, return previous
        //  for easy operation mode switching
        //operationHandling handleOperationMode(const operationHandling newOp);

        //- Return Well name
        const word& name() const {
            return name_;
        }

        //- Return model hosting dictionary
        const dictionary& dict() const {
            return wellDict_;
        }

        //- Return injection phase name
        const word& iPhase() const {
            return iPhase_;
        }

        //- Return well operation mode as enum
        const operationHandling& operation() const {
            return operation_;
        }

        //- Return Well total volume
        const scalar& totalVolume() const {
            return tV_;
        }

        //- Return cells set
        cellSet& wellSet() {
            return wellSet_;
        }

        //- Const access to cells IDs
        const labelList cellIDs() const {
            return wellSet_.toc();
        }

        //- Read perforations
        //void readPerforations();

        //- Const access to next timeStep to catch rate changes
        const scalar& nextTimeStep() const {
            return timeForDt_;
        }

        //- Correct well coeffs
        virtual void correct();

        //- Setup before correcting well coeffs
        virtual void preCorrect() = 0;

        //- Cleanup after correcting well coeffs
        virtual void postCorrect() = 0;

        //- Return target flowRate
        virtual const dimensionedScalar targetRate() const = 0;

};

class wellModel
{

protected:

    // Protected data

        //- Model Name
        word name_;

        //- Well dictionary entry
        dictionary wellProperties_;

        //- Const ref to mesh
        const fvMesh& mesh_;

        //- Const ref to pressure field
        const volScalarField& p_;

public:

    //- Runtime type information
    TypeName("wellModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        wellModel,
        dictionary,
        (
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh
        ),
        (name, wellProperties, mesh)
    );

    // Selectors

        //- Return a reference to the well control model
        static autoPtr<wellModel> New
        (
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh
        );


    // Constructors

        //- Construct from components
        wellModel
        (
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~wellModel() {}


    // Member Functions

        //- Return well model name
        const word& name() const {
            return name_;
        }

        //- Return Well dict entry
        const dictionary& dict() const {
            return wellProperties_;
        }

        // Correct well rates
        virtual void correct() = 0;

        virtual void operator()(word& w) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Macro for RTS-enabled class initiated from this class template
#define makeRelativePermeabilityModelBase(KType, nPhases)                      \
                                                                               \
defineNamedTemplateTypeNameAndDebug(base##nPhases####KType##Well, 0);          \
defineTemplateRunTimeSelectionTable(base##nPhases####KType##Well, dictionary); \

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "wellBaseAndModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
