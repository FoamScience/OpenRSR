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
#include "incompressibleFluid.H"
#include "relativePermeabilityModelBase.H"
#include "relativePermeabilityModelBasesFwd.H"
#include "capillaryPressureModelBase.H"
#include "capillaryPressureModelBasesFwd.H"
#include "wellBaseAndModel.H"
#include "wellBasesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
using namespace phaseModels;

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
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {

        Info<< "Time = " << runTime.timeName() << nl << endl;
        runTime++;

        // Solve saturation equation explicitly
        #include "alphaEqn.H"
        ////scalar minSw = gMin(phasew.alpha());
        ////if (minSw < 0)
        ////{
        ////    phasew.alpha() -= minSw;
        ////}
        ////phasew.alpha().correctBoundaryConditions();
        #include "updateSaturationFields.H"

        //// Solve pressure equation (implicit)
        #include "pEqn.H"
        #include "updatePressureFields.H"


        #include "setDeltaT.H"
        #include "updateWellFields.H"

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
