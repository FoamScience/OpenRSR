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
    Foam::wellModelBase

Description
    An abstract class for well models based on the assumption that
    well bottom hole pressure is a function of (main) pressure
    and some flowrate variable (Total, or phase-specific).

SourceFiles
    wellBase.C

\*---------------------------------------------------------------------------*/

#ifndef wellBase_H
#define wellBase_H

#include "phase.H"
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
#include "interpolationTable.H"
#include "wellCompileCalculations.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class wellBase Declaration
\*---------------------------------------------------------------------------*/

template<class KType,int nPhases> class wellModelBase;
class phase;

template<class KType, int nPhases>
class wellBase
    : public regIOobject
{

public:

    // Public data members

        //- Enumeration for handling operation modes
        enum operationHandling
        {
            PROD,             /*!< Production well */
            INJE,             /*!< Injection well */
        };

        //- A Hash table for handling imposed drives
        //  Currently reservers room to phase rates and BHP
        //  Todo: add room for phaseDissolvedInPhaseRatio
        HashTable<int> driveHandling;

protected:

    // Protected data

        //- Well Name
        word name_;

        //- Well configs dict
        const dictionary& wellDict_;

        //- Const-ref to mesh
        const fvMesh& mesh_;

        //- Const-ref to the model driving this well
        const wellModelBase<KType, nPhases>& corrector_;

        //- Well operation mode: injection - production
        operationHandling operation_;

        //- Well Orientation
        word orientation_;

        //- Injection phase name
        word iPhase_;

        //- Well radius
        dimensionedScalar radius_;

        //- Well skin
        scalar skin_;

        //- Well Perforations
        PtrList<topoSetSource> perfos_;

        //- Well topo set
        cellSet wellSet_;

        //- Phase names list
        wordList phases_;

        //- Well source coeffs as a function of pressure
        //  (Explicit and implicit parts)
        List<VectorN<scalar, 2*nPhases> > source_;

        //- Well total volume
        scalar  tV_;

        //- Well BHP
        dimensionedScalar bhp_;

        //- Target well rate
        dimensionedScalar tRate_;

        //- Imposed Drive series
        List<interpolationTable<scalar> > driveSeries_;

        //- TimeStep to catch imposed wellRate/BHP changes
        scalar timeForDt_;

    // Protected Member Functions

        //- Disallow copy construct
        wellBase(const wellBase&);

        //- Disallow default bitwise assignment
        void operator=(const wellBase&);

        //- Calculate total cell volume
        inline void cellsVolume()
        {
            tV_ = 0;
            forAll(wellSet_.toc(), ci)
            {
                tV_ += mesh_.V()[wellSet_.toc()[ci]];
            }
        }


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
            const dictionary& wellDict,
            const fvMesh& mesh,
            const wellModelBase<KType, nPhases>& corrector
        ),
        (name, wellDict, mesh, corrector)
	); 		

    // Constructors

        //- Construct from components
        wellBase
        (
            const word& name,
            const dictionary& wellDict,
            const fvMesh& mesh,
            const wellModelBase<KType, nPhases>& corrector
        );

    // Selectors

        //- Return a reference to the selected Kr model
        static autoPtr<wellBase<KType, nPhases> > New
        (
            const word& name,
            const dictionary& wellDict,
            const fvMesh& mesh,
            const wellModelBase<KType, nPhases>& corrector
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

        //- Return target flowRate
        const dimensionedScalar& targetRate() const
        {
            return tRate_;
        }

        //- Return skin value
        const scalar& skin() const
        {
            return skin_;
        }

        //- Return well radius
        scalar radius() const
        {
            return radius_.value();
        }

        //- Set well source coeffs
        void setSource(List<VectorN<scalar,4> >& src)
        {
            source_ = src;
        }

        //- Return Const-Ref to well Source Matrix
        const List<VectorN<scalar, 2*nPhases> >& source() const
        {
            return source_;
        }

        //- Check if a drive is active
        bool isActiveDrive(const word& drive) const
        {
            if (!driveHandling.found(drive))
                FatalErrorIn(__PRETTY_FUNCTION__)
                    << "Didn't expect drive: " << drive
                    << " at all for well: " << name_
                    << ".\n List of legitimate well drives include:"
                    << driveHandling.toc() << exit(FatalError);
            if (driveHandling[drive] > 0)
            {
                return 1;
            } else {
                return 0;
            }
        }

        //- Retrieve imposed drive value from series
        scalar driveAtTime(const word& activeDrive, scalar index) const
        {
            scalar imposedDrive = -1;
            if (isActiveDrive(activeDrive))
            {
                imposedDrive = driveSeries_[driveHandling[activeDrive]-1](index);
            }
            return imposedDrive;
        }

        //- Read perforations
        void readPerforations();

        //- Read imposed well rate/BHP
        void readImposedDrives();

        //- Const access to next timeStep to catch rate changes
        const scalar& nextTimeStep() const {
            return timeForDt_;
        }

        //- Figure out lower cell index
        label lowerCell() const {
            label minCell = cellIDs()[0];
            forAll(cellIDs(), ci)
            {
                scalar cg = (mesh_.C()[ci] && corrector_.g()).value();
                scalar base = (mesh_.C()[minCell] && corrector_.g()).value();
                minCell = (cg >= base) ? cellIDs()[ci] : minCell;
            }
            return minCell;
        }

        //- Figure out upper cell index
        label upperCell() const {
            label maxCell = cellIDs()[0];
            forAll(cellIDs(), ci)
            {
                scalar cg = (mesh_.C()[ci] && corrector_.g()).value();
                scalar base = (mesh_.C()[maxCell] && corrector_.g()).value();
                maxCell = (cg <= base) ? cellIDs()[ci] : maxCell;
            }
            return maxCell;
        }

        //- Correct well coeffs
        virtual void correct();

        //- Pre-Setup for well correction
        virtual void preCorrect() = 0;

        //- Post-Cleanup for well correction
        virtual void postCorrect() = 0;
};

template<class KType, int nPhases>
class wellModelBase
: public objectRegistry
{

protected:

    // Protected data

        //- Model Name
        const word& name_;

        //- Well dictionary entry
        const dictionary& wellProperties_;

        //- Const ref to mesh
        const fvMesh& mesh_;

        //- Const ref to pressure field
        const volScalarField& p_;

        //- Const ref to absolute permeability
        const KType& K_;

        //- Wells List
        PtrList<wellBase<KType, nPhases> > wells_;

        //- Sources matrices
        HashTable<fvScalarMatrix> source_;

        //- Ref to Gravitational field
        const uniformDimensionedVectorField& g_;

    // Potected Member Functions

        //- Read wells from list
        void readWells();

public:

    //- Runtime type information
    TypeName("wellBaseModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        wellModelBase,
        Dictionary,
        (
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh
        ),
        (name, wellProperties, mesh)
    );

    // Selectors

        //- Return a reference to the well control model
        static autoPtr<wellModelBase> New
        (
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh
        );


    // Constructors

        //- Construct from components
        wellModelBase
        (
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~wellModelBase() {}


    // Member Functions

        //- Return well model name
        const word& name() const {
            return name_;
        }

        //- Return Wells dict
        const dictionary& dict() const {
            return wellProperties_;
        }

        //- Return the object registry
        const objectRegistry& thisDb() const
        {
            return *this;
        }

        //- Return wells List
        const PtrList<wellBase<KType, nPhases> >& wells() const {
            return wells_;
        }

        //- Return phase well source matrix
        const fvScalarMatrix& source(const word& phaseName) const {
            return source_[phaseName];
        }

        //- Return a ref to model gracitational field
        const uniformDimensionedVectorField& g() const {
            return g_;
        }

        //- Correct well rates
        virtual void correct() = 0;

        //- Virtual Functor operator
        virtual void operator()(const word& wellName) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Macro for RTS-enabled class initiated from this class template
#define makeWellBase(KType, nPhases)                      \
                                                                               \
defineNamedTemplateTypeNameAndDebug(base##nPhases####KType##Well, 0);          \
defineTemplateRunTimeSelectionTable(base##nPhases####KType##Well, dictionary); \

#define makeWellModelBase(KType, nPhases)                      \
                                                                               \
defineNamedTemplateTypeNameAndDebug(base##nPhases####KType##WellModel, 0);          \
defineTemplateRunTimeSelectionTable(base##nPhases####KType##WellModel, Dictionary); \

#define makeWellModelDerived(Suffix, KType, MuRhoType)                      \
                                                                                \
defineNamedTemplateTypeNameAndDebug(derived##KType####MuRhoType##Suffix, 0);    \
addToRunTimeSelectionTable                                                      \
(                                                                               \
    base2##KType##WellModel,                                              \
    derived##KType####MuRhoType##Suffix,                                            \
    Dictionary                                                                  \
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "wellBaseAndModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
