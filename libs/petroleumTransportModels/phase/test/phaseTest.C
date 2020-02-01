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
    Test the base class "phase"

Developers
    F. Mohamed Elwardi

\*---------------------------------------------------------------------------*/

#define CATCH_CONFIG_MAIN
#include "catch.H"

#include "fvCFD.H"
#include "phase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

SCENARIO("Phase Object Creation with Default configuration", "[Virtual]")
{
    GIVEN("Time, mesh, a child of 'phase' class and a valid transportProperties dictionary")
    {
        #include "createTestTimeAndMesh.H"

        // A child phase object
        class childPhase : public phase
        {
            public:
                childPhase(
                    const word& name,
                    const dictionary& transportProperties,
                    const fvMesh& mesh
                ) : phase(name, transportProperties, mesh) { }
                virtual void correct(){
                    phiPtr_->correctBoundaryConditions();
                    U_ = fvc::reconstruct(phiPtr_());
                    U_.correctBoundaryConditions();
                }
                bool writeData(Ostream& os) const { notImplemented("bool writeData(Ostream&) const"); }
        };

        word phaseName = "water";

        // The transportProperties dict
        dictionary transportProperties;
        transportProperties.add(phaseName, dictionary());

        WHEN("Calling correct() on a stack-allocated phase object")
        {
            // Needs the presence of '0/water.U' dictionary
            childPhase waterPhase(phaseName, transportProperties, mesh);

            surfaceScalarField& phi = waterPhase.phi();
            volVectorField& U = waterPhase.U();

            // Emitate a parabolic U profile
            forAll(U.internalField(), ci)
            {
                U.internalField()[ci] = vector(mesh.nCells() - ci, 0, 0);
            }
            U.correctBoundaryConditions();

            // Needs to set div(water.phi, water.U) in system/fvSchemes
            waterPhase.correct();

            THEN("Velocity must be constructed properly from phase flux using Weller's approximation")
            {
                // Particularly test if fvc::reconstruct has changed!
                // fvc::reconstruct may be dissipative near some BCs
                surfaceScalarField newPhi = linearInterpolate(U) & mesh.Sf();
                REQUIRE( phi.internalField()==newPhi.internalField() );
            }
        }

        WHEN("Allocating a child phase object on the stack with default configuration")
        {
            // Needs the presence of '0/water.U' dictionary
            childPhase waterPhase(phaseName, transportProperties, mesh);

            THEN("Phase saturation must be equal to 1")
            {
                const volScalarField& Sw = waterPhase.alpha();
    
                // Require that all internal cell values are 1s
                REQUIRE( Sw.internalField() == scalarList(mesh.nCells(), 1.0));

                // Require that all (non-empty) boundary values are 1s
                const wordList& BCtypes = mesh.boundaryMesh().types();
                forAll(BCtypes, patchi)
                {
                    const scalarList& bSw = Sw.boundaryField()[patchi];
                    if (BCtypes[patchi] != "empty")
                    REQUIRE( bSw == scalarList(bSw.size(), 1.0) );
                }

                // Print (byte) size of phase object
                //Info << int(sizeof(waterPhase)) << endl;
            }
        }

    }
}

// ************************************************************************* //
