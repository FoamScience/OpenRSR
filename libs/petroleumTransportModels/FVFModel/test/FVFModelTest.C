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
    Test the base class "FVFModel"

Developers
    F. Mohamed Elwardi

\*---------------------------------------------------------------------------*/

#define CATCH_CONFIG_MAIN
#include "catch.H"

#include "fvCFD.H"
#include "FVFModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// A child FVFModel
// A local class can't have static members, so we put this class here

namespace Foam
{
    class childFVFModel : public FVFModel
    {
        public:
            TypeName("childFVFModel");
            childFVFModel(
                const word& name,
                const dictionary& phaseDict,
                const volScalarField& p
            ): FVFModel(name, phaseDict, p){};
            ~childFVFModel(){}
            void correct() {
                forAll(rFVF_.internalField(), ci){
                    rFVF_.internalField()[ci] = ci;
                }
                rFVF_.correctBoundaryConditions();
            }
    };
    defineTypeNameAndDebug(childFVFModel, 0);
    addToRunTimeSelectionTable(FVFModel, childFVFModel, dictionary);
}

using namespace Foam;


SCENARIO("FVF Model Selection with Default configuration", "[Virtual]")
{
    GIVEN("Pressure field and RunTime-Selectable child of 'FVFModel' class")
    {
        FatalError.throwExceptions();
        #include "createTestTimeAndMesh.H"

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

        // The phase dictionary, which has a 'FVFModel' entry
        dictionary dict = dictionary::null;
        dict.add(word("FVFModel"), word(""));

        WHEN("Passing invalid model identifier to the selector")
        {
            // Try to instantiate a model from an "InvalidModelName" class
            dict.set("FVFModel", "InvalidModelName");

            THEN("An exception is raised")
            {
                REQUIRE_THROWS(FVFModel::New("invalidModel", dict, p));
            }
        }

        WHEN("Calling correct() on a valid FVF model object")
        {
            // Try to instantiate a model from an "InvalidModelName" class
            dict.set("FVFModel", "childFVFModel");
            autoPtr<FVFModel> fvf = FVFModel::New("fvf", dict, p);
            fvf->correct();
            THEN("1|FVF Boundary conditions must reflect zeroGradient situation")
            {
                const volScalarField& rfvf = fvf->rFVF();
                forAll(mesh.boundary(), pi)
                {
                    const word& patchType = mesh.boundary()[pi].type();
                    if (patchType != "empty")
                    forAll(mesh.boundary()[pi], fi)
                    {
                        const label& bcell = mesh.boundaryMesh()[pi].faceCells()[fi];
                        CHECK(rfvf.internalField()[bcell] == rfvf.boundaryField()[pi][fi]);
                    }
                }
            }
        }

        WHEN("Calling correct() on a valid FVF model object")
        {
            // Try to instantiate a model from an "InvalidModelName" class
            dict.set("FVFModel", "childFVFModel");
            autoPtr<FVFModel> fvf = FVFModel::New("fvf", dict, p);
            fvf->correct();
            THEN("d(1|FVF)dP Boundary conditions must reflect zeroGradient situation")
            {
                const volScalarField& drfvfdp = fvf->rFVF();
                forAll(mesh.boundary(), pi)
                {
                    const word& patchType = mesh.boundary()[pi].type();
                    if (patchType != "empty")
                    forAll(mesh.boundary()[pi], fi)
                    {
                        const label& bcell = mesh.boundaryMesh()[pi].faceCells()[fi];
                        CHECK(drfvfdp.internalField()[bcell] == drfvfdp.boundaryField()[pi][fi]);
                    }
                }
            }
        }
    }
}

// ************************************************************************* //
