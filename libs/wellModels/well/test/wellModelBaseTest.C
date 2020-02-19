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
    wellModelBaseTest

Description
    Test the instantiated base classes for different well models.

Developers
    F. Mohamed Elwardi

\*---------------------------------------------------------------------------*/

#define CATCH_CONFIG_MAIN
#include "catch.H"

#include "fvCFD.H"
#include "wellBaseAndModel.H"
#include "parameterTypes.H"
#include "wellBasesFwd.H"
#include "incompressibleFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
using namespace phaseModels;

namespace Foam
{
    // Add a wellModelDerived class
    #include "testChildClasses/derived2PhasesIsoWellModel.H"
    // Add a childWell class
    #include "testChildClasses/derived2PhasesIsoWell.H"
}


SCENARIO("Construction of well objects")
{
    GIVEN("A valid mesh, and a wellProperties file containing a list of well entries")
    {
        #include "createTestTimeAndMesh.H"
        #include "testChildClasses/setupTwoPhaseCase.H"

        WHEN("An iso-2-phases well model object is created")
        {
            autoPtr<base2IsoWellModel>    wModel =
                base2IsoWellModel::New("derivedWellModel", wellProperties, mesh);

            PtrList<entry> wellsEntries = wellProperties.lookup("wells");
            CHECK(wModel->wells().size() == wellsEntries.size());
            CHECK (wModel->wells()[0].name() == word("aWell"));

            THEN("Perforation intervals are read correctly")
            {
                labelList wellCells(4);
                wellCells[0] = 0;
                wellCells[1] = 1;
                wellCells[2] = 5;
                wellCells[3] = 6;
                REQUIRE(wModel->wells()[0].cellIDs() == wellCells);
            }

            THEN("Imposed phase rates must be correctly read")
            {
                // Check that BHP drive is inactive
                CHECK(wModel->wells()[0].driveAtTime("BHP", 1000) == -1);
                REQUIRE(wModel->wells()[0].driveAtTime(water.name()+".rate", 1) == 0.6);
            }
        }

        WHEN("Pass a wellModel child to wellBase's constructor, and run wellBase.correct()")
        {
            autoPtr<base2IsoWellModel>    wModel =
                base2IsoWellModel::New("derivedWellModel", wellProperties, mesh);
            THEN("wellBase Members must be updated properly")
            {
                wModel->wells()[0].correct();
                REQUIRE(wModel->wells()[0].equivalentRadius() == 15);
            }
        }
    }
}

// ************************************************************************* //
