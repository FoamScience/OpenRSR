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
    relativePermeabilityModelBaseTest

Description
    Test the instantiated base classes base2PhasesKrModel and
    base3PhasesKrModel

Developers
    F. Mohamed Elwardi

\*---------------------------------------------------------------------------*/

#define CATCH_CONFIG_MAIN
#include "catch.H"

#include "fvCFD.H"
#include "relativePermeabilityModelBase.H"
#include "relativePermeabilityModelBasesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

SCENARIO("Kr-related compile time calculations")
{
    GIVEN("An int N, representing number of phases for Kr model")
    {
        WHEN("N is a Compile-Time property)")
        {
            const int N = 2;
            THEN("Number of Kr field must be correctly computed at compiled-time.")
            {
                CHECK(numberOfKrFields(2) == 2);
                CHECK(numberOfKrFields(3) == 6);
                // Logic: This is a compile-time test!
                // The size of the c-like "array" is 1
                // which has to be known at compile-time
                // So the hole condition evaluation takes place at compile time
                // if the number of Kr fields is valid, the branch throwing exception
                // is never evaluated, and the array has a size of 1 (a const)
                const char array[numberOfKrFields(N) == 2 ? 1 : throw std::exception() ];
            }
        }
    }
}

namespace Foam
{
// A child relativePermeabilityModelBase class
// providing an implementation of correct()
class childKrModel : public relativePermeabilityModelBase<2>
{
public:
    TypeName("childKrModel")
    childKrModel(
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
    ): relativePermeabilityModelBase<2>(name, transportProperties, mesh){}
    virtual ~childKrModel(){}
    virtual void correct(){}
    virtual const wordList phaseNames() const {return wordList();};
    virtual const word canonicalPhase() const {return word("");  };
    label phaseIndex(const word& phaseName) const
    {return 0;};
};

#include "addToRunTimeSelectionTable.H"

defineTypeNameAndDebug(childKrModel, 0);
addToRunTimeSelectionTable(base2PhasesKrModel, childKrModel, dictionary);

}


SCENARIO("Base relative permeability class template")
{
    GIVEN("A valid mesh, and a transport properties file")
    {
        #include "createTestTimeAndMesh.H"

        word phaseName = "oil";

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

        WHEN("A 2-Phases Kr Model is instantiated with krTable of value uniform (1 0)")
        {
            autoPtr<base2PhasesKrModel> krModel = base2PhasesKrModel::New("krModel", transportProperties, mesh);
            THEN("Access to kr tables works as expected")
            {
                // Canonical phase Kr
                volScalarField kr0 = krModel->kr(0);
                std::vector<double> ones(mesh.nCells(), 1);
                std::vector<double> kr0Vec(mesh.nCells(), 0);
                forAll(kr0Vec, ci)
                {
                    kr0Vec[ci] = kr0.internalField()[ci];
                }
                REQUIRE(kr0Vec == ones);
            }
        }
    }
}

// ************************************************************************* //
