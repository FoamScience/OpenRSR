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
    phaseTest

Description
    Tests for the incomressible phase object

Developers
    F. Mohamed Elwardi

\*---------------------------------------------------------------------------*/

#define CATCH_CONFIG_MAIN
#include "catch.H"

#include "fvCFD.H"
#include "incompressibleFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
using namespace phaseModels;

SCENARIO("Incompressible Phase object interface")
{
    GIVEN("Time, mesh, and a valid phase dictionary (providing rho, mu)")
    {
        #include "createTestTimeAndMesh.H"

        word phaseName = "water";

        // The transportProperties dict (with rho and mu)
        dictionary transportProperties;

        // A sub dict with rho - mu as dimensionedScalars
        dictionary phaseDict;
        phaseDict.add("rho", dimensionedScalar(phaseName, dimMass/dimVolume, 1000));
        phaseDict.add("mu", dimensionedScalar(phaseName, dimMass/dimLength/dimTime, 1e-3));

        // Add the subdict to the parent transportDict
        transportProperties.add(phaseName, phaseDict);

        WHEN("correct() is called on an on-the-stack incompressibleFluid")
        {
            // Needs the presence of '0/water.U' dictionary
            incompressibleFluid waterPhase(phaseName, transportProperties, mesh);
            waterPhase.correct();

            THEN("Retreived Rho value and dimensions must be consistent")
            {
                REQUIRE( waterPhase.rho().value() == 1000.0 );
                REQUIRE( waterPhase.rho().dimensions() == dimMass/dimVolume );
            }
        }

        WHEN("correct() is called on an on-the-stack incompressibleFluid")
        {
            // Needs the presence of '0/water.U' dictionary
            incompressibleFluid waterPhase(phaseName, transportProperties, mesh);
            waterPhase.correct();

            THEN("Retreived Mu value and dimensions must be consistent")
            {
                REQUIRE( waterPhase.mu().value() == 1e-3 );
                REQUIRE( waterPhase.mu().dimensions() == dimMass/dimLength/dimTime );

                // Size of phase, and incompressibleFluid
                //Info << int(sizeof(phase)) << " " << int(sizeof(waterPhase)) << endl;
            }
        }

    }
}

// ************************************************************************* //
