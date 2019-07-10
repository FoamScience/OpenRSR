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
    impesFoam

Description
    Transient solver for incompressible two-phase flow (Darcy's law) in porous media
    using the IMPES method (IMplicit Pressure Explicit Saturation).
    Permeability is isotropic (K == volScalarField)

Developers
    F. Mohamed Elwardi

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "harmonic.H"
#include "fvBlockMatrix.H"
#include "incompressibleFluid.H"
#include "capillaryPressureModel.H"
#include "relativePermeabilityModel.H"
#include "wellModel.H"
//#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "delta.H"
    #include "createTimeControls.H"
    #include "readGravitationalAcceleration.H"

    #include "createFields.H"
    #include "createSaturationsFields.H"
    #include "createPressureFields.H"
    #include "createWellFields.H"

    #include "readTimeControls.H"
    #include "readSaturationBounds.H"
    #include "initContinuityErrs.H"
    #include "initCNVErrors.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        #include "readBlockSolverControls.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        fvBlockMatrix<vector2> resEqn(Swp);
        #include "SEqn.H"
        #include "pEqn.H"
        resEqn.solve();
        //Info<< "Solver performance: " << Swp.name() 
        //    << " Initial Residual " << prf.initialResidual() 
        //    << ", Final Residual " << prf.finalResidual() 
        //    << ", NoIters " << prf.nIterations() 
        //    << endl;

        resEqn.retrieveSolution(0, phasew.alpha().internalField());
        resEqn.retrieveSolution(1, p.internalField());
        #include "boundSaturation.H"

        phasew.alpha().correctBoundaryConditions();
        p.correctBoundaryConditions();

        #include "CNVErrors.H"
        #include "updateSaturationFields.H"
        #include "updatePressureFields.H"
        #include "updateWellFields.H"

        Info<< "Saturation w " << " Min = " << gMin(phasew.alpha()) 
            << " Max = " << gMax(phasew.alpha()) <<  endl;

        runTime.write();
        #include "setDeltaT.H"
        runTime++;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
