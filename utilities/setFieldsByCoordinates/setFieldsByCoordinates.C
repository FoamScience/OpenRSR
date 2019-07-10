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
    setFieldByCoordinates

Description
    A tool to take away the burden of properly "IDing" OpenFOAM mesh
    cells to match MUFITS/ECLIPSE ones! It's expected to generate a table
    of cell center coordinate along with field(s) data from Paraview.
    Supports only one volScalarField at a time for now.

    Drived with a Dict where list of fields and their data are specified.
    psiData ((xCoord yCoord zCoord psi));

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
    Info<< "Reading setFieldsByCoordinatesDict\n" << endl;
    IOdictionary dict
    (
        IOobject
        (
            "setFieldsByCoordinatesDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Fields name
    const wordList psiNames = dict.lookup("fields");

    // Treat each specified field seperately
    forAll(psiNames, fieldi){
        // Name of current field
        const word psiName = psiNames[fieldi];
        // List of data
        const List<vector4> psiData = dict.lookup(psiName+"Data");

        // Field in question
        volScalarField psi
        (
            IOobject
            (
                psiName,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("",dimless,VGREAT),
            "zeroGradient"
        );

        // Fatal error if sizes don't mach
        if (psiData.size() != psi.internalField().size()) 
            WarningIn("int main(int argc, char *argv[])") 
                << psiName << "has " << psi.internalField().size() << "cells but" 
                << psiData.size() << "enountered in " 
                << psiName+"Data entry." << endl;


        // Set field values
        forAll(psiData, celli){

            vector cellCoords(psiData[celli][0],psiData[celli][1],psiData[celli][2]);
            label cellID = mesh.findCell(cellCoords);
            if (cellID == -1) 
                FatalErrorIn("int main(int argc, char *argv[])") 
                    << "No cell in mesh for point " << cellCoords << endl;
            psi[cellID] = psiData[celli][3];
        }

        // Check for unvisited cells
        bool unvisited = true;
        while(unvisited)
        {
            unvisited = false;
            forAll(psi.internalField(), celli)
            {
                if (psi[celli] > Foam::sqrt(VGREAT))
                {
                    // Cell missing in original data
                    // Volume averaging
                    labelList neighbors = mesh.cellCells()[celli];
                    scalar psiV = 0;
                    forAll(neighbors, nei)
                    {
                        psiV += mesh.V()[neighbors[nei]]*psi[neighbors[nei]];
                    }
                    psi[celli] = psiV/mesh.V()[celli];
                    unvisited = true;
                }
            }
        }

        // write the field
        Info<< "Writing field " << psiName << endl;
        psi.write();

    }
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
