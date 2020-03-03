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

#include "verticalIsoWell.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace twoPhaseWells
{

defineTypeNameAndDebug(verticalIsoWell, 0);
addToRunTimeSelectionTable(base2IsoWell, verticalIsoWell, dictionary);


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

verticalIsoWell::verticalIsoWell
(
    const word& name,
    const dictionary& wellProperties,
    const fvMesh& mesh,
    const base2IsoWellModel& corrector
)
:
    base2IsoWell(name, wellProperties, mesh, corrector)
{
}

// * * * * * * * * * * * * * * Public Member Methods * * * * * * * * * * * * * //

bool verticalIsoWell::writeData(Ostream&) const
{
}

void verticalIsoWell::preCorrect()
{
}

void verticalIsoWell::postCorrect()
{
}

};  // End namespace twoPhaseWellModels

}  // End namespace Foam


