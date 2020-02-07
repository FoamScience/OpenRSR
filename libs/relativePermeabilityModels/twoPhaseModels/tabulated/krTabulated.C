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

\*---------------------------------------------------------------------------*/

#include "krTabulated.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{
    defineTypeNameAndDebug(krTabulated, 0);
    addToRunTimeSelectionTable
    (
        base2PhasesKrModel,
        krTabulated,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePermeabilityModels::krTabulated::krTabulated
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
    :
    base2PhasesKrModel(name, transportProperties, mesh),
    coeffsDict_(transportProperties.subDict(typeName+"Coeffs")),
    phaseNames_
    (
        coeffsDict_.found("phases")
        ? coeffsDict_.lookup("phases")
        : mesh.objectRegistry::names<phase>()
    ),
    canonicalPhase_(coeffsDict_.lookup("canonicalPhase")),
    S_(mesh.objectRegistry::lookupObject<volScalarField>(canonicalPhase_+".alpha")),
    krSeries_(coeffsDict_)
{
}

// * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * * //

void Foam::relativePermeabilityModels::krTabulated::correct()
{
    forAll(mesh_.C(), celli)
    {
        // Calculate Krs
        kr_.internalField()[celli][0] = krSeries_(S_.internalField()[celli])[0];
        kr_.internalField()[celli][1] = krSeries_(S_.internalField()[celli])[1];

        // Calculate Kr Derivatives irt S_
        dkrdS_.internalField()[celli][0] = krSeries_(S_.internalField()[celli])[2];
        dkrdS_.internalField()[celli][1] = krSeries_(S_.internalField()[celli])[3];
    }

    //forAll(kr1_.boundaryField(), patchi)
    //{
    //    forAll(kr1_.boundaryField()[patchi], facei)
    //    {
    //        scalar Satf =  S_.boundaryField()[patchi][facei];
    //        kr1_.boundaryField()[patchi][facei] = krSeries_(Satf)[0];
    //        kr2_.boundaryField()[patchi][facei] = krSeries_(Satf)[1];
    //        dkr1dS_.boundaryField()[patchi][facei] = krSeries_(Satf)[2];
    //        dkr2dS_.boundaryField()[patchi][facei] = krSeries_(Satf)[3];
    //    }
    //}

    // It's unbelievable how important these lines are !!!!
    kr_.correctBoundaryConditions();
    dkrdS_.correctBoundaryConditions();

}

// ************************************************************************* //
