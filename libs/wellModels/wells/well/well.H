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
    A namespace to hold different kinds of well objects and their positioning
    methods.

Class
    Foam::well

Description
    An abstract base class for well objects. Keeps track of all well objects
    in a static list and registers each well object to mesh.objectRegistry for
    individual processing.

SourceFiles
    well.C

\*---------------------------------------------------------------------------*/

#ifndef well_H
#define well_H

#include "cellSet.H"
#include "HCmixture.H"
#include "volFields.H"
#include "regIOobject.H"
#include "blackoilPhase.H"
#include "topoSetSource.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class well Declaration
\*---------------------------------------------------------------------------*/

class well
    : public regIOobject
{

public:

    // Public data members
        
        //- Global well objects list
        static PtrList<well> wellsList;

        //- Global well names list
        static wordList wellsNames;

        //- Enueration for handling well flow modes
        enum rateHandling
        {
            LRATE,             /*!< Constant Total flowRate */
            //ORATE,             /*!< Constant oil flowRate */
            //WRATE,             /*!< Constant water flowRate */
            BHP,             /*!< Constant BHP */
        };

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

        //- Well dictionary entry
        dictionary wellDict_;

        //- Ref to mesh
        const fvMesh& mesh_;

        //- Ref to phase objects
        const HCmixture<blackoilPhase>& mixture_;

        //- Well Type: vertical - horizontal
        word type_;

        //- Well control mode: flowRate - BHP
        rateHandling mode_;

        //- Well operation mode: injection - production
        operationHandling operation_;

        //- Injection phase
        word iPhase_;

        //- well radius
        dimensionedScalar radius_;

        //- Perforations
        PtrList<topoSetSource> perfos_;

        //- well topo set
        cellSet wellSet_;

       //- well total volume
       scalar  tV_;

       struct flowRateModeVars
       {
           //- Well (target) Bottom Hole flowRate
           interpolationTable<scalar> flowRateSeries;

           // Constructors

               //- Default Constructor
               flowRateModeVars(){}

               //- Construct from Components
               flowRateModeVars(const dictionary& wellDict);

       } flowRateData_;

       //- Well water & oil source coefficients
       scalarList ewSource_,eoSource_;
       scalarList iwSource_,ioSource_;

       //- Are perforation intervals already read?
       bool perfosNotRead_;

       //- TimeStep to catch well rate changes
       scalar timeForDt_;

       //- Allow switching operation modes
       bool allowSwitchingModes_;

       //- Is this well switched already to BHP-operated
       bool switchedToBHP_;

        //- Matrix for the WJ.p term
        fvScalarMatrix WJw_;


public:

    //- Runtime type information
    TypeName("well");


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        well,
        dictionary,
        (
            const word& name,
            const dictionary& wellDict,
            const fvMesh& mesh,
            const HCmixture<blackoilPhase>& mixture
        ),
        (name, wellDict, mesh, mixture)
    );


    // Selectors

        //- Return a reference to the selected permeability model
        static autoPtr<well> New
        (
            const word& name,
            const dictionary& wellDict,
            const fvMesh& mesh,
            const HCmixture<blackoilPhase>& mixture
        );

        //- Return a reference to the selected permeability model
        static void readWells
        (
            const dictionary& wellsProperties,
            const fvMesh& mesh,
            const HCmixture<blackoilPhase>& mixture
        );

    // Constructors

        //- Construct from components
        well
        (
            const word& name,
            const dictionary& wellDict,
            const fvMesh& mesh,
            const HCmixture<blackoilPhase>& mixture
        );


    //- Destructor
    virtual ~well() {}


    // Member Functions

        //- Return operationHandling as word
        word operationHandlingToWord(const operationHandling& op) const;

        //- Return rateHandling as word
        word rateHandlingToWord(const rateHandling& rate) const;

        //- Return operationHandling as enumeration
        operationHandling wordToOperationHandling(const word& op) const;

        //- Return rateHandling as enumeration
        rateHandling wordToRateHandling(const word& rate) const;

        //- Handle operation mode from enum, return previous
        //  for easy operation mode switching
        operationHandling handleOperationMode(const operationHandling newOp);

        //- Handle flow mode from enum, return previous
        //  for easy rate control switching
        rateHandling handleRateMode(const rateHandling& newRate);

        //- Return Well name
        const word& name() const {
            return name_;
        }

        //- Return Well type
        const word& wellType() const {
            return type_;
        }

        //- Return injection phase name
        const word& iPhase() const {
            return iPhase_;
        }

        //- Return controlMode as enum
        const rateHandling& mode() const {
            return mode_;
        }

        //- Return well operation mode as enum
        const operationHandling& operation() const {
            return operation_;
        }

        //- Return Well dict entry
        const dictionary& dict() const {
            return wellDict_;
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
        const dimensionedScalar targetRate() const {
            if (mode_ == well::LRATE or switchedToBHP_)
            {
                scalar t = mesh_.time().timeOutputValue();
                return dimensionedScalar
                    (name_+"_flowRate",dimVolume/dimTime, flowRateData_.flowRateSeries(t));
            } else {
                FatalErrorIn("const dimensionedScalar Foam::well::flowRate() const;")
                    << "Well " << name_ << " is not operating in flowRate mode."
                    << exit(FatalError);
            }
        }

        //- Const access to explicit oil cell sources
        const scalarList eoSource() const {
            return eoSource_;
        }

        //- Const access to explicit water cell sources
        const scalarList ewSource() const {
            return ewSource_;
        }

        //- Const access to implicit oil cell sources
        const scalarList ioSource() const {
            return ioSource_;
        }

        //- Const access to implicit water cell sources
        const scalarList iwSource() const {
            return iwSource_;
        }

        //- Const access to perforations status
        bool perfosNotRead() {
            return perfosNotRead_;
        }

        void readPerforations();

        //- Const access to next timeStep to catch rate changes
        const scalar& nextTimeStep() const {
            return timeForDt_;
        }

        //- Const access to next timeStep to catch rate changes
        const bool& allowSwitchingModes() const {
            return allowSwitchingModes_;
        }

        //- Switch rate-imposed and BHP-imposed modes in a well 
        virtual void switchWellMode() = 0;

        ////- Estimate equivalent radius for the well
        //virtual void estimateEquivRadius() = 0;

        ////- Calculate well index 
        //virtual void calculateWellIndex() = 0;

        //- Calculate well sources
        //virtual void sourceWell() = 0;

        //- Calculate well sources
        virtual void correct() = 0;

        //- Return WJ matrix for switched wells
        const fvScalarMatrix& WJmatrix() const
        {
            return WJw_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
