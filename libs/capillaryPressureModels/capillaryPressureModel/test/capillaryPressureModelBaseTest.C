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
    capillaryPressureModelBaseTest

Description
    Test the instantiated base classes base2PhasesPcModel and
    base3PhasesPcModel

Developers
    F. Mohamed Elwardi

\*---------------------------------------------------------------------------*/

#define CATCH_CONFIG_MAIN
#include "catch.H"

#include "fvCFD.H"
#include "capillaryPressureModelBase.H"
#include "capillaryPressureModelBasesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

namespace Foam
{
// A child capillaryPressureModelBase class
// providing an implementation of correct()
class childPcModel : public capillaryPressureModelBase<2>
{
public:
    TypeName("childPcModel")
    childPcModel(
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
    ): capillaryPressureModelBase<2>(name, transportProperties, mesh){}
    virtual ~childPcModel(){}
    virtual void correct(){}
    virtual const wordList phaseNames() const {return wordList();};
    virtual const word canonicalPhase() const {return word("");  };
};

#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(childPcModel, 0);
addToRunTimeSelectionTable(base2PhasesPcModel, childPcModel, dictionary);

}


SCENARIO("Base capillarity class template")
{
    GIVEN("A valid mesh, and a transport properties file")
    {
        #include "createTestTimeAndMesh.H"

        // The transportProperties dict
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

        WHEN("A 2-Phases Pc Model is instantiated with pcTable of value uniform 0.5")
        {
            autoPtr<base2PhasesPcModel> pcModel = base2PhasesPcModel::New("pcModel", transportProperties, mesh);
            THEN("Access to pc fields works as expected")
            {
                // Canonical phase Pc
                auto pc0 = pcModel->pcTable();
                std::vector<double> ones(mesh.nCells(), 0.5);
                std::vector<double> pc0Vec(mesh.nCells(), 0);
                forAll(pc0Vec, ci)
                {
                    pc0Vec[ci] = pc0.internalField()[ci];
                }
                REQUIRE(pc0Vec == ones);
            }
        }
    }
}

// ************************************************************************* //
