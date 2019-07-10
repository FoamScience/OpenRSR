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
    setAnisoPermeability

Description
    A tool to constuct the anisotropic Permeability field from 
    its three components at cell centers.

Developers
    F.M. Elwardi

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createTimeControls.H"
    #include "readTimeControls.H"


    // A driving dictionary 
    Info<< "Reading setAnisoPermeabilityDict\n" << endl;
    IOdictionary dict
    (
        IOobject
        (
            "setAnisoPermeabilityDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // The field
    volTensorField psi
    (
        IOobject
        (
            "K",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("",dimArea,
            tensor(VSMALL,VSMALL,VSMALL,
                   VSMALL,VSMALL,VSMALL,
                   VSMALL,VSMALL,VSMALL)),
        "zeroGradient"
    );

    // The data list In the format
    // (X Y Z Kxx Kyy Kzz)
    const List<vector6> KData = dict.lookup("KData");

    // Set the field
    forAll(KData, celli){

        vector cellCoords(KData[celli][0],KData[celli][1],KData[celli][2]);
        label cellID = mesh.findCell(cellCoords);
        if (cellID == -1) 
            FatalErrorIn("int main(int argc, char *argv[])") 
                << "No cell in mesh for point " << cellCoords << endl;
            psi[cellID][0] = KData[celli][3];
            psi[cellID][4] = KData[celli][4];
            psi[cellID][8] = KData[celli][5];
    }

    // Check for still default cells
    forAll(psi.internalField(), celli)
    {
        if (psi[celli][0] <= VSMALL 
                or psi[celli][4] <= VSMALL
                or psi[celli][8] <= SMALL)
        {
            Info << nl << "Didn't visit cell " << celli;
        }
    }

    // write the field
    Info<< nl << nl << "Writing field K" << endl;
    psi.write();

    return 0;
}
