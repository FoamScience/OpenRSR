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

#include "Peaceman.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blackoilWellModels
{
defineTypeNameAndDebug(Peaceman, 0);
addToRunTimeSelectionTable
(
    wellModel,
    Peaceman,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blackoilWellModels::Peaceman::Peaceman
(
    const word& name,
    const dictionary& wellsProperties,
    const fvMesh& mesh,
    const HCmixture<blackoilPhase>& mixture,
    const relativePermeabilityModel& krModel
)
:
    wellModel(name, wellsProperties, mesh, mixture, krModel)
{
}

void Foam::blackoilWellModels::Peaceman::correct()
{
    // Update Diag/source coeffs
    forAll(wells_, welli)
    {
        wells_[welli].correct();    
        const labelList& cells = wells_[welli].wellSet().toc();

        forAll(cells, celli)
        {
            eoSource_[cells[celli]] = wells_[welli].eoSource()[celli];
            ewSource_[cells[celli]] = wells_[welli].ewSource()[celli];
            ioSource_[cells[celli]] = wells_[welli].ioSource()[celli];
            iwSource_[cells[celli]] = wells_[welli].iwSource()[celli];
        }

        const scalar& wellTimeForDt = wells_[welli].nextTimeStep();
        if (timeForDt_ > wellTimeForDt)
        {
            timeForDt_ = wellTimeForDt;
        }
    }

}

// ************************************************************************* //
