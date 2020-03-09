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

Application
    peacemanTest

Description
    Test Two-Phase Peaceman well model features

Developers
    F. Mohamed Elwardi

\*---------------------------------------------------------------------------*/

#define CATCH_CONFIG_MAIN
#include "catch.H"

#include "fvCFD.H"
#include "phase.H"
#include "incompressibleFluid.H"
#include "wellBaseAndModel.H"
#include "wellBasesFwd.H"
#include "peaceman.H"
#include "relativePermeabilityModelBase.H"
#include "relativePermeabilityModelBasesFwd.H"
#include "capillaryPressureModelBase.H"
#include "capillaryPressureModelBasesFwd.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
using namespace phaseModels;

SCENARIO("Instantiation of two-Phase Peaceman well model object")
{
    GIVEN("Time, A valid mesh, Kr and Pc models")
    {
        #include "createTestTimeAndMesh.H"
        // transportProperties dict
        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        // The wellProperties dict
        IOdictionary wellProperties
        (
            IOobject
            (
                "wellProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        // The pressure field
        volScalarField p
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("p", dimPressure, 1e7)
        );
        // The permeability field
        volScalarField K
        (
            IOobject
            (
                "K",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("K", dimArea, 1e-12)
        );

        // Gravitation field
        #include "readGravitationalAcceleration.H"

        // Two phases to work with
        incompressibleFluid oil(word("oil"), transportProperties, mesh);
        incompressibleFluid water(word("water"), transportProperties, mesh);

        // Setup Kr and Pc model
        autoPtr<base2PhasesKrModel> krModel = 
            base2PhasesKrModel::New("krModel", transportProperties, mesh);
        autoPtr<base2PhasesPcModel> pcModel = 
            base2PhasesPcModel::New("pcModel", transportProperties, mesh);

        WHEN("A 2-cells well is selected from a Peaceman model")
        {
            autoPtr<base2IsoWellModel> wModel =
                base2IsoWellModel::New("wModel", wellProperties, mesh);

            const base2IsoWell& aWell = 
                wModel->objectRegistry::lookupObject<base2IsoWell>("aWell");
            CHECK(aWell.cellIDs().size() == 2);
            THEN("Lower cell is correctly selected")
            {
                REQUIRE(aWell.lowerCell() == 0);
            }

            THEN("Upper cell is correctly selected")
            {
                REQUIRE(aWell.upperCell() == 1);
            }
        }

        WHEN("Total canonical-phase flow rate is supplied for a multi-cell well")
        {
            autoPtr<base2IsoWellModel> wModel =
                base2IsoWellModel::New("wModel", wellProperties, mesh);

            const base2IsoWell& aWell = 
                wModel->objectRegistry::lookupObject<base2IsoWell>("aWell");
            CHECK(aWell.cellIDs().size() == 2);

            THEN("The total well rate reported is correct.")
            {
                const double totalRateFromFile = aWell.driveAtTime("water.rate", 1);
                double calculatedTotalRate = 0;

                runTime.setTime(1.0, 1);
                CHECK(aWell.isActiveDrive("water.rate"));
                krModel->correct();
                wModel->correct();

                const fvScalarMatrix& wellMatrix = wModel->source("water");
                if (wModel->debug())
                {
                    Info << wellMatrix.source() << nl
                         << "Is asymmetric? " << wellMatrix.asymmetric() << nl
                         << "Is diagonal? " << wellMatrix.asymmetric() << nl
                         << "Dimensions: " << wellMatrix.dimensions() << nl
                         << "Lower: " << wellMatrix.hasLower() << nl
                         << "Upper: " << wellMatrix.hasUpper() << nl;
                }

                //calculatedTotalRate = sum(wellMatrix.source());
                forAll(wellMatrix.source(), ci)
                {
                    calculatedTotalRate += wellMatrix.source()[ci]*mesh.V()[ci];
                }
                REQUIRE_THAT(totalRateFromFile,
                        Catch::WithinAbs(-aWell.operationSign()*calculatedTotalRate, 0.001)
                );
            }
        }
    }
}

// ************************************************************************* //
