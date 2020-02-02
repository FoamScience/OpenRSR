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
    FVFModelTest

Description
    Test the tabulated FVF class

Developers
    F. Mohamed Elwardi

\*---------------------------------------------------------------------------*/

#define CATCH_CONFIG_MAIN
#include "catch.H"

#include "fvCFD.H"
#include "blackoilPhase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

SCENARIO("Black oil Phase Interface integrity")
{
    GIVEN("Time, mesh, and a valid phase dictionary (providing standard rho, FVF model)")
    {
        #include "createTestTimeAndMesh.H"

        word phaseName = "oil";

        // The transportProperties dict (with rho and mu)
        dictionary transportProperties;

        // Some levels of sub-dicts
        dictionary phaseDict = dictionary::null;
        dictionary FVFDict = dictionary::null;

        phaseDict.add(word("rhoSc"), dimensionedScalar("rhoSc",dimMass/dimVolume,0.7e3));
        phaseDict.add(word("FVFModel"), word("tabulated"));

        dictionary fvfData = dictionary::null;
        FVFDict.add("fileName", "\"FVF.dat\"");
        FVFDict.add("outOfBounds", "warn");

        // Add the subdicts to their parents 
        phaseDict.add(word("FVFData"), FVFDict);
        transportProperties.add(phaseName, phaseDict);


        WHEN("correct() is called while mesh has a record of a pressure field (p)")
        {
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
                dimensionedScalar("p", dimPressure, 1.8e7)
            );

            // Construct a phase object
            blackoilPhase oilPhase(phaseName, transportProperties, mesh);

            THEN("Density must be correctly updated according to FVF and p.")
            {
                oilPhase.correct();

                std::vector<scalar> rhoCalc(mesh.nCells(), 0);
                std::vector<scalar> rhoManual(mesh.nCells(), 540.91964);
                forAll(rhoCalc, ci)
                {
                    rhoCalc[ci] = oilPhase.rho().internalField()[ci];
                }

                REQUIRE_THAT( rhoCalc, Catch::Matchers::Approx(rhoManual) );
            }
        }
    }
}

// ************************************************************************* //
