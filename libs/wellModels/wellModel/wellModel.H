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
    Foam::blackoilWellModels

Description
    A namespace for well models: classes that actually return matrix
    coefficients for well sources in blackoil equations

Class
    Foam::wellModel

Description
    An abstract base class for wellModels. Constructs matrix coeffs for
    wells in BlackOil simulations. Works only for two-phase simulations.

SourceFiles
    wellModel.C

\*---------------------------------------------------------------------------*/

#ifndef wellModel_H
#define wellModel_H

#include "well.H"
#include "fvMatrix.H"
#include "HCmixture.H"
#include "blackoilPhase.H"
#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class wellModel Declaration
\*---------------------------------------------------------------------------*/

class wellModel
{

protected:

    // Protected data
    
        //- Model Name
        word name_;

        //- Well dictionary entry
        dictionary wellsProperties_;

        //- Const ref to mesh
        const fvMesh& mesh_;

        //- Const ref to pressure field
        const volScalarField& p_;

        //- Ref to phase objects
        const HCmixture<blackoilPhase>& mixture_;

        //- Ref to relative permeability model
        const relativePermeabilityModel& krModel_;

        //- Ref to wells list
        PtrList<well>& wells_;

        //- Fields for semi-implicit well sources 
        //  The well source is then for oil: 
        //  "fvm::Sp(iosource, p) + eoSource()"
        volScalarField eoSource_, ioSource_;
        volScalarField ewSource_, iwSource_;

        //- Next time any rate-operated well changes its imposed rate
        scalar timeForDt_;

        //- Matrix for the WJ.p term
        fvScalarMatrix WJw_;


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
            const dictionary& wellsProperties,
            const fvMesh& mesh,
            const HCmixture<blackoilPhase>& mixture,
            const relativePermeabilityModel& krModel
        ),
        (name, wellsProperties, mesh, mixture, krModel)
    );


    // Selectors

        //- Return a reference to the well control model
        static autoPtr<wellModel> New
        (
            const word& name,
            const dictionary& wellsProperties,
            const fvMesh& mesh,
            const HCmixture<blackoilPhase>& mixture,
            const relativePermeabilityModel& krModel
        );


    // Constructors

        //- Construct from components
        wellModel
        (
            const word& name,
            const dictionary& wellsProperties,
            const fvMesh& mesh,
            const HCmixture<blackoilPhase>& mixture,
            const relativePermeabilityModel& krModel
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
            return wellsProperties_;
        }

        //- Return next timeStep to catch imposed rate changes
        scalar& nextTimeStep() {
            return timeForDt_;
        }

        volScalarField& ExplicitOilSource() {
            return eoSource_;
        }
        volScalarField& ImplicitOilSource() {
            return ioSource_;
        }
        volScalarField& ExplicitWaterSource() {
            return ewSource_;
        }
        volScalarField& ImplicitWaterSource() {
            return iwSource_;
        }

        // Correct well rates
        virtual void correct() = 0;

        //- Specialized treatment of wells
        void updateWellProperties();

        //- Return WJ matrix
        fvScalarMatrix& WJmatrix()
        {
            forAll(wells_, welli)
            {
                const well& thisWell = wells_[welli];
                WJw_ += thisWell.WJmatrix();
            }
            return WJw_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
