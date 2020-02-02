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
#include "FVFModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

SCENARIO("Tabulated FVF Model with before-and-after bubble point data")
{
    GIVEN("Pressure field and FVF data file (FVF.dat for example)")
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

        // The phase dictionary, which has a 'FVFModel' entry and FVFData subdict
        dictionary dict = dictionary::null;
        dict.add(word("FVFModel"), word("tabulated"));

        dictionary fvfData = dictionary::null;
        fvfData.add("fileName", "\"FVF.dat\"");
        fvfData.add("outOfBounds", "warn");
        dict.add("FVFData", fvfData);

        WHEN("Supplied valid P values")
        {
            autoPtr<FVFModel> fvf = FVFModel::New("fvf", dict, p);
            std::vector<scalar> pCheckList = {
                7.50e4, 1.00e5, 3.00e5, 4.00e5, 8.00e5, 1.50e6, 7.00e6,
                1.00e7, 1.80e7, 2.30e7
            };
            std::vector<scalar> rFVFManual = {
                0.94641579, 0.944130711, 0.925850035, 0.916709697, 0.885988783,
                0.86582205, 0.791810943, 0.769912425, 0.7798852, 0.786176065
            };
            std::vector<scalar> drFVFdPManual = {
                -9.13394e-8, -8.94851e-8, -7.46509e-8, -6.72338e-8, -3.50993e-8,
                -1.93656e-8, -8.70729e-9, 1.22407e-9, 1.43534e-9, 1.46173e-9
            };

            forAll(p.internalField(), ci)
            {
                p.internalField()[ci] = pCheckList[ci];
            }

            // Perform model calculations
            fvf->correct();

            std::vector<scalar> rFVFCalculated(pCheckList.size());
            forAll(rFVFCalculated, ci)
            {
                rFVFCalculated[ci] = fvf->rFVF().internalField()[ci];
            }
            std::vector<scalar> drFVFdPCalculated(pCheckList.size());
            forAll(drFVFdPCalculated, ci)
            {
                drFVFdPCalculated[ci] = fvf->drFVFdP().internalField()[ci];
            }

            THEN("Linear interpolation must be fairly accurate")
            {
                REQUIRE_THAT(rFVFCalculated, Catch::Matchers::Approx(rFVFManual));
                REQUIRE_THAT(drFVFdPCalculated, Catch::Matchers::Approx(drFVFdPManual));
            }
        }
    }
}

// ************************************************************************* //
