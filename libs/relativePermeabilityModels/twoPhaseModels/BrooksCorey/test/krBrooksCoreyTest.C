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
    krBrooksCoreyTest

Description
    Test Modified Brooks Corey Relative Permeability class

Developers
    F. Mohamed Elwardi

\*---------------------------------------------------------------------------*/

#define CATCH_CONFIG_MAIN
#include "catch.H"

#include "fvCFD.H"
#include "krBrooksCorey.H"
#include "incompressibleFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
using namespace phaseModels;

SCENARIO("Relative Permeability Curves in Water-Oil systems")
{
    GIVEN("A mesh, a couple of phases and model coeffs in transportProperties")
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
        // Two Black-Oil phases
        incompressibleFluid oil(word("oil"), transportProperties, mesh);
        incompressibleFluid water(word("water"), transportProperties, mesh);

        // The Kr Model
        autoPtr<base2PhasesKrModel> krModel = base2PhasesKrModel::New("krModel", transportProperties, mesh);

        THEN("correct() must calculate Kr values correctly")
        {
            const scalar Swi = 0.2;
            const scalar Snr = 0.1;
            forAll(mesh.C(), ci)
            {
                water.alpha().internalField()[ci] = Swi+ci*(1-Snr-Swi)/(mesh.nCells()-1);
            }
            
            krModel->correct();
            volScalarField krw = krModel->kr(0);
            volScalarField krn = krModel->kr(1);

            const std::vector<scalar> KrwManual = {
                0, 0.012345679, 0.049382716, 0.111111111,
                0.197530864, 0.308641975, 0.444444444, 0.604938272,
                0.790123457, 1
            };
            const std::vector<scalar> KrnManual = {
                0.9, 0.711111111, 0.544444444, 0.4, 0.277777778,
                0.177777778, 0.1, 0.044444444, 0.011111111, 0
            };
            std::vector<scalar> KrCalculated(mesh.nCells());
            forAll(KrCalculated, ci)
            {
                KrCalculated[ci] = krw.internalField()[ci];
            }
            REQUIRE_THAT(KrCalculated, Catch::Matchers::Approx(KrwManual).margin(1e-6));
            forAll(krn, ci)
            {
                KrCalculated[ci] = krn.internalField()[ci];
            }
            REQUIRE_THAT(KrCalculated, Catch::Matchers::Approx(KrnManual).margin(1e-6));
        }

        THEN("correct() must calculate Kr derivatives correctly")
        {
            const scalar Swi = 0.2;
            const scalar Snr = 0.1;
            forAll(mesh.C(), ci)
            {
                water.alpha().internalField()[ci] = Swi+ci*(1-Snr-Swi)/(mesh.nCells()-1);
            }
            
            krModel->correct();
            volScalarField dkrw = krModel->dkrdS(0);
            volScalarField dkrn = krModel->dkrdS(1);

            const std::vector<scalar> dKrwManual = {
                0, 0.317460317, 0.634920635, 0.952380952, 1.26984127,
                1.587301587, 1.904761905, 2.222222222, 2.53968254, 2.857142857
            };
            const std::vector<scalar> dKrnManual = {
                -2.571428571, -2.285714286, -2, -1.714285714, -1.428571429,
                -1.142857143, -0.857142857, -0.571428571, -0.285714286, 0
            };
            std::vector<scalar> dKrCalculated(mesh.nCells());
            forAll(dKrCalculated, ci)
            {
                dKrCalculated[ci] = dkrw.internalField()[ci];
            }
            REQUIRE_THAT(dKrCalculated, Catch::Matchers::Approx(dKrwManual).margin(1e-6));
            forAll(dKrCalculated, ci)
            {
                dKrCalculated[ci] = dkrn.internalField()[ci];
            }
            REQUIRE_THAT(dKrCalculated, Catch::Matchers::Approx(dKrnManual).margin(1e-6));
        }

    }
}

// ************************************************************************* //
