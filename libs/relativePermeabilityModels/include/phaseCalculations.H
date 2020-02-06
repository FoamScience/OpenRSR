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

Description
    (Allegedly) Compile-Time calculations for how many kr fields should 
    be there as a function of the number of phases involved.

\*---------------------------------------------------------------------------*/

#ifndef phaseCalculations_H
#define phaseCalculations_H

#include "VectorN.H"
#include "volMesh.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// CompileTime factorials
constexpr int cmpFactorial(int i)
{
    return i == 0 ? 1 : (i*cmpFactorial(i-1));
}

// Calculate number of Kr fields needed
// Number of Kr fields when "nPhases" phases are present
// That is, 2*numberOfPhaseInterfaces
// (Unique Permutations)
// Note: Forced compile-time execution
static constexpr int numberOfKrFields(int nPhases)
{
    return cmpFactorial(nPhases)/cmpFactorial(nPhases-2);
}

// An alias for Kr table field
template<int nPhases>
using volKrField = GeometricField< VectorN<scalar,nPhases>, fvPatchField, volMesh>;

// An alias for an undimensioned entry in Kr table field
template<int nPhases>
using krVector = VectorN<scalar, numberOfKrFields(nPhases)>;

// An alias for a dimensioned entry in Kr table field
template<int nPhases>
using dimensionedKrVector = dimensioned<VectorN<scalar, numberOfKrFields(nPhases)> >;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
