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

#include "HCmixture.H"

#ifndef HCmixture_C
#define HCmixture_C

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class FluidModel>
Foam::HCmixture<FluidModel>::HCmixture
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
:  
    regIOobject(
        IOobject
        (
            name,
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    HCmixtureDict_(transportProperties),
    mesh_(mesh),
    phaseNames_(HCmixtureDict_.lookup("phases")),
    wettingStatus_(phaseNames_.size()),
    phaseList_(phaseNames_.size()),
    pbModel_
    (
        bubblePointModel::New
        (
            "pbModel",
            transportProperties,
            mesh.lookupObject<volScalarField>("p")
        )
    )
{
    // Get wettability configuration
    forAll(wettingStatus_, phasei){
        wettingStatus_[phasei] = 
            HCmixtureDict_.subDict(phaseNames_[phasei]).lookupOrDefault<word>("wettingStatus", "wetting");
    }

    // Create pointers to phases
    forAll(phaseList_, phasei){
        phaseList_[phasei] = 
            FluidModel::New
            (
                phaseNames_[phasei],
                HCmixtureDict_,
                mesh_
            );
    }
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class FluidModel>
Foam::autoPtr<Foam::HCmixture<FluidModel> > Foam::HCmixture<FluidModel>::New
(
    const word& name,
    const dictionary& transportProperties,
    const fvMesh& mesh
)
{
    return autoPtr<HCmixture<FluidModel> >
    (
        new HCmixture(name, transportProperties, mesh)
    );
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class FluidModel>
inline FluidModel& Foam::HCmixture<FluidModel>::operator[](label i)
{
    return phaseList_[i]();
}

#endif
