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
    (Allegedly) Compile-Time calculations for well classes

\*---------------------------------------------------------------------------*/

#include "krPcCompileCalculations.H"
#include <vector>

namespace Foam
{

//- Template Function to find binary combinations
template<class T>
constexpr List<List<T> > findBinaryCombinations(const List<T>& a)
{
    const int N = a.size();
    List<List<T> > permsList(int(cmpFactorial(N)/2));

    std::vector<bool> v(N);
    std::fill(v.begin(), v.begin() + 2, true);

    int ii = 0;
    do {
        for (int i = 0; i < N; ++i) {
            List<bool> vList(N, true);
            forAll(vList, vi)
            {
                vList[vi] = v[vi];
            }
            permsList[ii] = subset(vList, a);
        }
        ii++;
    } while (std::prev_permutation(v.begin(), v.end()));
   return permsList;
}

}
