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
    wellBaseTest

Description
    Test the instantiated base classes for different types of wells.

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

struct FunctorTestDerived : public wellModel
{
    TypeName("FunctorTestDerived");
    FunctorTestDerived(
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh
            ): wellModel(name, wellProperties, mesh) {}
    virtual void correct(){}
    virtual void operator()(const word& wellName) const
    {
        if (foundObject<base2IsoWell>(wellName))
        {
            const base2IsoWell& well = lookupObject<base2IsoWell>(wellName);
            Info << "Got const well ref: " << well.cellIDs() << endl;
        }
    }
};
#include "addToRunTimeSelectionTable.H"
defineTypeNameAndDebug(FunctorTestDerived, 0);
addToRunTimeSelectionTable
(
    wellModel,
    FunctorTestDerived,
    dictionary
);

// A child wellBase class
// providing an implementation of correct()
class childWell : public wellBase<Iso,2>
{
public:
    TypeName("childWell")
    childWell(
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh,
            const wellModel& corrector
    ): wellBase<Iso,2>(name, wellProperties, mesh, corrector){}
    virtual ~childWell(){}
    virtual bool writeData(Ostream&) const {}
    virtual void preCorrect(){ Info << "Executing preCorrect"  << endl; }
    virtual void postCorrect(){ Info << "Executing postCorrect"  << endl; }
};
defineTypeNameAndDebug(childWell, 0);
addToRunTimeSelectionTable(base2IsoWell, childWell, dictionary);

}


SCENARIO("Integration between wellBase class templates and wellModel class")
{
    GIVEN("A valid mesh, and a wellProperties file")
    {
        #include "createTestTimeAndMesh.H"

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
        // The transportProperties dict (with rho and mu)
        dictionary transportProperties;

        // A sub dict with rho - mu as dimensionedScalars
        dictionary waterDict;
        waterDict.add("rho", dimensionedScalar(word("water"), dimMass/dimVolume, 1000));
        waterDict.add("mu", dimensionedScalar(word("water"), dimMass/dimLength/dimTime, 1e-3));
        // A sub dict with rho - mu as dimensionedScalars
        dictionary oilDict;
        oilDict.add("rho", dimensionedScalar(word("oil"), dimMass/dimVolume, 700));
        oilDict.add("mu", dimensionedScalar(word("oil"), dimMass/dimLength/dimTime, 1e-5));

        // Add the subdict to the parent transportDict
        transportProperties.add(word("water"), waterDict);
        // Add the subdict to the parent transportDict
        transportProperties.add(word("oil"), waterDict);


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
            dimensionedScalar("p", dimPressure, 0.0)
        );

        // Two phases to work with
        incompressibleFluid oil(word("oil"), transportProperties, mesh);
        incompressibleFluid water(word("water"), transportProperties, mesh);

        WHEN("Pass a wellModel child to wellBase's constructor, and run wellBase.correct()")
        {
            autoPtr<wellModel>    wModel = 
                wellModel::New("derivedWellModel", wellProperties, mesh);
            autoPtr<base2IsoWell> aWell  = 
                base2IsoWell::New( "aWell", wellProperties, mesh, wModel());
            THEN("wellBase Members must be updated properly")
            {
                CHECK(aWell->iPhase() == word("noPhase"));
                aWell->correct();
                //REQUIRE(aWell->iPhase() == word("noPhaseTestFunctor"));
            }
        }

        WHEN("A Well object is constructed")
        {
            autoPtr<wellModel>    wModel = 
                wellModel::New("derivedWellModel", wellProperties, mesh);
            autoPtr<base2IsoWell> aWell  = 
                base2IsoWell::New( "aWell", wellProperties, mesh, wModel());
            THEN("Perforation intervals are read correctly")
            {
                CHECK(wellProperties.subDict("aWell").found("perforations"));
                labelList wellCells(4);
                wellCells[0] = 0;
                wellCells[1] = 1;
                wellCells[2] = 5;
                wellCells[3] = 6;
                REQUIRE(aWell->cellIDs() == wellCells);
            }
        }

        WHEN("A 2-phases well object is constructed")
        {
            autoPtr<wellModel>    wModel = 
                wellModel::New("derivedWellModel", wellProperties, mesh);
            autoPtr<base2IsoWell> aWell  = 
                base2IsoWell::New( "aWell", wellProperties, mesh, wModel());
            THEN("Imposed phase rates must be correctly read")
            {
                // Check that BHP drive is inactive
                CHECK(aWell->driveAtTime("BHP", 1000) == -1);
                REQUIRE(aWell->driveAtTime(water.name()+".rate", 1) == 0.6);
            }
        }
    }
}

// ************************************************************************* //
