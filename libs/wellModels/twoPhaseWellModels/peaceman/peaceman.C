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

#include "peaceman.H"
#include "fvmSup.H"


namespace Foam
{

//namespace twoPhaseWellModels
//{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class KType>
peaceman<KType>::peaceman
(
    const word& name,
    const dictionary& wellProperties,
    const fvMesh& mesh
)
:
    wellModelBase<KType,2>(name, wellProperties, mesh),
    krModel_
    (
        mesh.objectRegistry::lookupObject<relativePermeabilityModelBase<2> >
        ("krModel")
    ),
    pcModel_
    (
        mesh.objectRegistry::lookupObject<capillaryPressureModelBase<2> >
        ("pcModel")
    ),
    canonicalMu_
    (
        mesh.objectRegistry::lookupObject<volScalarField>(
            krModel_.canonicalPhase()+".mu"
        )
    ),
    nonCanonicalMu_
    (
        mesh.objectRegistry::lookupObject<volScalarField>(
            (krModel_.phaseNames()[1] == krModel_.canonicalPhase() )
            ? krModel_.phaseNames()[0]+".mu"
            : krModel_.phaseNames()[1]+".mu"
        )
    ),
    canonicalRho_
    (
        mesh.objectRegistry::lookupObject<volScalarField>(
            krModel_.canonicalPhase()+".rho"
        )
    ),
    nonCanonicalRho_
    (
        mesh.objectRegistry::lookupObject<volScalarField>(
            (krModel_.phaseNames()[1] == krModel_.canonicalPhase() )
            ? krModel_.phaseNames()[0]+".rho"
            : krModel_.phaseNames()[1]+".rho"
        )
    )
{
}

// * * * * * * * * * * * * * * Private Member Methods * * * * * * * * * * * * * //

template<class KType>
scalarList peaceman<KType>::calculateFractionalFlow
(
    const wellBase<KType, 2>& well
)
{
    //  Where Qi = -(K kri A/mu)*grad(p)
    // This ignores the effect of Pc on fractional flow for now
    // Logic: grad(pc) should be much less than grad(p)
    // Needs more investigation though
    scalarList fracFlow(well.cellIDs().size());
    forAll(fracFlow, ci)
    {
        label cellID = well.cellIDs()[ci];
        fracFlow[ci] = 1.0
            /(1.0 + krModel_.kr(1, cellID)*canonicalMu_[cellID]
                /(krModel_.kr(0, cellID)*nonCanonicalMu_[cellID]+SMALL)
            );
    }
    return fracFlow;
}

template<class KType>
List<vector> peaceman<KType>::estimateCellSizes
(
    const wellBase<KType, 2>& well
)
{
    List<vector> h(well.cellIDs().size());
    // First estimate cell sizes
    forAll(h, ci)
    {
        label cellID = well.cellIDs()[ci];
        const labelList& cellEdges = this->mesh_.cellEdges()[cellID];
        scalarList dEdgesAlongZ(cellEdges.size());
        scalarList dEdgesAlongX(cellEdges.size());
        scalarList dEdgesAlongY(cellEdges.size());
        forAll(cellEdges, edgei){
            dEdgesAlongX[edgei] = mag(this->mesh_.edges()[cellEdges[edgei]].vec(this->mesh_.points()).x());
            dEdgesAlongY[edgei] = mag(this->mesh_.edges()[cellEdges[edgei]].vec(this->mesh_.points()).y());
            dEdgesAlongZ[edgei] = mag(this->mesh_.edges()[cellEdges[edgei]].vec(this->mesh_.points()).z());
        }
        h[ci][0] = gMax(dEdgesAlongX);
        h[ci][1] = gMax(dEdgesAlongY);
        h[ci][2] = gMax(dEdgesAlongZ);
    }
    return h;
}

template<class KType>
scalarList peaceman<KType>::estimateEquivRadius
(
    const wellBase<KType, 2>& well
)
{
    notImplemented(__PRETTY_FUNCTION__);
}

template<>
scalarList peaceman<Iso>::estimateEquivRadius
(
    const wellBase<Iso, 2>& well
)
{
    // Isotropic media ---> Re depends only on geometry
    scalarList re(well.cellIDs().size());

    if (well.dict().found("equivalentRadius"))
    {
        // If re list is provided, use it
        // (For optimization and history matching)
        well.dict().lookup("equivalentRadius") >> re;
    } else {
        // If not provided, use a standard estimation
        List<vector> h = estimateCellSizes(well);
        forAll(re, ci)
        {
            re[ci] = 0.14*Foam::sqrt(Foam::pow(h[ci][0], 2) + Foam::pow(h[ci][1], 2));
        }
    }
    return re;
}

template<class KType>
scalarList peaceman<KType>::calculateWellPI
(
    const Foam::wellBase<KType, 2>& well
)
{
    notImplemented(__PRETTY_FUNCTION__);
}


template<>
scalarList peaceman<Iso>::calculateWellPI
(
    const Foam::wellBase<Iso, 2>& well
)
{
    List<vector> h = estimateCellSizes(well);
    scalarList re = estimateEquivRadius(well);
    scalarList pi(re.size());

    forAll(pi, ci)
    {
        label cellID = well.cellIDs()[ci];
        pi[ci] = 2*mathematicalConstant::pi*K_[cellID]*h[ci][2]
            /(log(re[ci]/well.radius()) + well.skin());
    }
    return pi;
}

// * * * * * * * * * * * * * * Public Member Methods * * * * * * * * * * * * * //

template<class KType>
void peaceman<KType>::correct()
{
    // Boundary conditions are completely ignored for now
    dimensionedScalar dzero("canIm", dimless/dimTime/dimPressure, 0.0);
    volScalarField canIm
    (
        IOobject
        (
            "canIm",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dzero,
        "zeroGradient"
    );
    volScalarField canEx
    (
        IOobject
        (
            "canEx",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dzero,
        "zeroGradient"
    );
    volScalarField ncaIm
    (
        IOobject
        (
            "ncaIm",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dzero,
        "zeroGradient"
    );
    volScalarField ncaEx
    (
        IOobject
        (
            "ncaEx",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dzero,
        "zeroGradient"
    );
    
    forAll(this->wells_, wi)
    {
        this->wells_[wi].correct();
        forAll(this->wells_[wi].cellIDs(), ci)
        {
            label cellID  = this->wells_[wi].cellIDs()[ci];
            canIm[cellID] = this->wells_[wi].source()[ci][0];
            canEx[cellID] = this->wells_[wi].source()[ci][1];
            ncaIm[cellID] = this->wells_[wi].source()[ci][2];
            ncaEx[cellID] = this->wells_[wi].source()[ci][3];
        }
    }

    this->source_[0] = fvm::Sp(canIm, this->p_) + fvm::Su(canEx, this->p_);
    this->source_[1] = fvm::Sp(ncaIm, this->p_) + fvm::Su(ncaEx, this->p_);
}

// * * * * * * * * * * * * * * * Public Operators * * * * * * * * * * * * * * * //

template<class KType>
void peaceman<KType>::operator()(const word& wellName) const
{
//    // Get const-ref to requested well
//    const wellBase<KType, 2>& well = 
//        objectRegistry::lookupObject<wellBase<KType,2> >(wellName);
//    // Update Well Productivity Index
//    scalarList pi = estimateWellPI(well);
//
//    // Setup well source matrix
//    // Thats finding coeffs where:
//    // Qi = ai*BHP + bi*P + ci
//    // for each phase
//    List<VectorN<scalar,4> > src(pi.size());
//    scalar timeIndex = 
//        wellModelBase<KType,2>::mesh_.time().timeOutpuValue();
//    word nonCanonicalPhase =
//            (krModel_.phaseNames()[1] == krModel_.canonicalPhase() )
//            ? krModel_.phaseNames()[0] 
//            : krModel_.phaseNames()[1];
//
//    if(well.isActiveDrive(krModel_.canonicalPhase()+".rate"))
//    {
//        // Case 1: canonical Q is constrained
//        // Estimate Canonical Fractional Flow
//        scalarList fQ = calculateFractionalFlow(well);
//        forAll(src, ci)
//        {
//            label cellID = well.cellIDs()[ci];
//            // P and BHP coeffs are null
//            src[ci][0] = 0;
//            src[ci][2] = 0;
//            // Free term is constrained constant for each phase
//            // Use fractional flow to express non canonical rate
//            src[ci][1] = well.driveAtTime(krModel_.canonicalPhase()+".rate", timeIndex);
//            src[ci][3] = src[ci][1]*(1-fQ[cellID])/fQ[cellID];
//        }
//    } else if (well.isActiveDrive(nonCanonicalPhase+".rate"))
//    {
//        // Case 2: non-canonical Q is constrained
//        notImplemented("nonCanonical Q constrainte not avilable")
//
//    } else if (well.isActiveDrive("totalRate"))
//    {
//        // Case 3: Total Q is constrained
//        notImplemented("TotalRate constrainte not avilable")
//    } else if (well.isActiveDrive("BHP"))
//    {
//        // Case 4: BHP is constrained
//        notImplemented("BHP constrainte not avilable")
//    }
//
//    // Update well source matrix
//    well.setSource(src);
}

//};  // End namespace twoPhaseWellModels

}  // End namespace Foam


