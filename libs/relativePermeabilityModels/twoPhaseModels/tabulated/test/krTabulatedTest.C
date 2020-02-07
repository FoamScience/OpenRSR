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
    krTabulatedTest

Description
    Test Tabulated Kr data entry

Developers
    F. Mohamed Elwardi

\*---------------------------------------------------------------------------*/

#define CATCH_CONFIG_MAIN
#include "catch.H"

#include "fvCFD.H"
#include "krTabulated.H"
#include "incompressibleFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
using namespace phaseModels;

SCENARIO("Relative Permeability Curves in Water-Oil systems")
{
    GIVEN("A mesh, a couple of phases and model configs in transportProperties")
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
        // Two incompressible phases
        incompressibleFluid oil(word("oil"), transportProperties, mesh);
        incompressibleFluid water(word("water"), transportProperties, mesh);

        // The Kr Model
        autoPtr<base2PhasesKrModel> krModel = base2PhasesKrModel::New("krModel", transportProperties, mesh);

        // The Kr table should reflect 3rd degree MBC model

        THEN("correct() must calculate Kr values correctly")
        {
            water.alpha().internalField()[0] = 0.3;

            krModel->correct();
            volScalarField krw = krModel->kr(0);
            volScalarField krn = krModel->kr(1);

			REQUIRE_THAT(krw.internalField()[0], Catch::Matchers::WithinAbs(0.002915452, 6e-3));
			REQUIRE_THAT(krn.internalField()[0], Catch::Matchers::WithinAbs(0.566763848, 6e-3));
        }

        THEN("correct() must calculate Kr derivatives correctly near the end-points")
        {
            water.alpha().internalField()[0] = 0.88;

            krModel->correct();
            volScalarField dkrw = krModel->dkrdS(0);
            volScalarField dkrn = krModel->dkrdS(1);

			REQUIRE_THAT(dkrw.internalField()[0], Catch::Matchers::WithinAbs(4.044314869, 2e-2));
			REQUIRE_THAT(dkrn.internalField()[0], Catch::Matchers::WithinAbs(-0.003148688, 2e-2));
        }

    }
}

// ************************************************************************* //
