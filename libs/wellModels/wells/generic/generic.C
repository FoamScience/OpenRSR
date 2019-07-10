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

#include "generic.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class PermeabilityType>
Foam::wells::generic<PermeabilityType>::generic
(
    const word& name,
    const dictionary& wellsProperties,
    const fvMesh& mesh,
    const HCmixture<blackoilPhase>& mixture
)
    :
    well(name, wellsProperties, mesh, mixture),
    kr1_(mesh.lookupObject<volScalarField>(mixture.phase("wetting").name()+".kr")),
    kr2_(mesh.lookupObject<volScalarField>(mixture.phase("nonWetting").name()+".kr")),
    mu1_(mesh.lookupObject<volScalarField>(mixture.phase("wetting").name()+".mu")),
    mu2_(mesh.lookupObject<volScalarField>(mixture.phase("nonWetting").name()+".mu")),
    rho1_(mesh.lookupObject<volScalarField>(mixture.phase("wetting").name()+".rho")),
    rho2_(mesh.lookupObject<volScalarField>(mixture.phase("nonWetting").name()+".rho")),
    waterFrac_(mesh.lookupObject<volScalarField>("Fw")),
    g_
    (
        mesh.lookupObject<uniformDimensionedVectorField>("g")
    ),
    BHPData_(),
    TfName_(wellDict_.lookupOrDefault<word>("TfName","Tof")),
    muName_(wellDict_.lookupOrDefault<word>("muName","oil.mu")),
    p_(mesh.lookupObject<volScalarField>("p")),
    Tf_(mesh.lookupObject<surfaceScalarField>(TfName_)),
    mu_(mesh.lookupObject<volScalarField>(muName_)),
    meshSubset_
    (
        IOobject
        (
            name_+".meshSubset",
            mesh_.time().constant(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    cellmap_(),
    facemap_(),
    submesh_(mesh)
{
    // Conditionally initialize well control variables
    if (mode_ == well::LRATE)
    {
        flowRateData_ = flowRateModeVars(wellDict_);
        Info<< "Well " << name_
            << " is operating in " << mode_ << "mode."
            << endl;
    } else if (mode_ == well::BHP)
    {
        BHPData_ = BHPModeVars(name_,wellDict_,mesh_);
        Info<< "Well " << name_
            << " is operating in " << mode_ << "mode."
            << endl;
    } else {
        FatalErrorIn("Foam::well::well(const word& name, const dictionary& wellDict, const fvMesh& mesh)")
            << "Well control mode was not recognized; available conrol modes: 2(flowRate BHP);" 
            << exit(FatalError);
    }
    ewSource_.setSize(wellSet_.toc().size());

    readPerforations();
    meshSubset_.setLargeCellSubset(wellSet_),
    cellmap_ = meshSubset_.cellMap();
    facemap_ = meshSubset_.faceMap();
    const fvMesh& submesh_(meshSubset_.subMesh());
    tP_.reset
    ( 
    new volScalarField
        (
            IOobject
            (
                name_+".pTest",
                submesh_.time().timeName(),
                submesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            submesh_,
            dimensionedScalar("pRef",dimPressure,0.0)
        )
    );
    tQ_.reset
    ( 
    new volScalarField
        ( 
            IOobject
            (
                name_+".qTest",
                submesh_.time().timeName(),
                submesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            submesh_,
            dimensionedScalar("qRef",dimless/dimTime,0.0),
            "zeroGradient"
        )
    );
    tTf_.reset
    ( 
    new surfaceScalarField
        (
            IOobject
            (
                name_+".Tf",
                submesh_.time().timeName(),
                submesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            submesh_,
            dimensionedScalar("TfRef",dimPressure,0.0),
            calculatedFvPatchScalarField::typeName
        )
    );
    if (tTf_().internalField().size() < 1)
    {
        FatalErrorIn("Foam::wells::generic<PermeabilityType>::generic\
                        (   \
                            const word& name, \
                            const dictionary& wellsProperties,\
                            const fvMesh& mesh,\
                            const HCmixture<blackoilPhase>& mixture\
                        )")
            << "Well orientation: generic doesn't work with single-cell wells"
            << " because you can't solve diffusivity equation on one cell!"
            << exit(FatalError);
    }
    forAll(cellmap_, celli)
    {
        tP_()[celli] = p_[cellmap_[celli]];
    }
    forAll(facemap_, facei)
    {
        tTf_()[facei] = Tf_[facemap_[facei]];
    }

    // fix flux bug with foam-extend
    submesh_.schemesDict().setFluxRequired(tP_().name());
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

template<class PermeabilityType>
void Foam::wells::generic<PermeabilityType>::calculateWellIndex()
{
    // Map global fields to well
    forAll(cellmap_, celli)
    {
        tP_()[celli] = p_[cellmap_[celli]];
    }
    forAll(facemap_, facei)
    {
        tTf_()[facei] = Tf_[facemap_[facei]];
    }

    // Compressibility
    dimensionedScalar comp("comp", dimless/dimPressure, 1e-14);

    // Solve diffusivity eqn in well mesh
    fvScalarMatrix wellMatrix
    (
        comp*fvm::ddt(tP_())+fvm::laplacian(-tTf_(), tP_())
    );
    wellMatrix.solve();
    tQ_() = fvc::div(wellMatrix.flux());

    Info << tQ_() << tP_() << nl;
    Info<< "######################"""" " << nl;

    // Calculate well PI at each well cell
    forAll(wellSet_, celli)
    {
        BHPData_.pi[celli] =
            tQ_()[celli]/(p_[cellmap_[celli]]-tP_()[celli]);
    }
    // Debug code
    if(1)
    {
        Info<< "Well " << name_ << " :" << nl
            << "PI:" << tab << BHPData_.pi << endl;
    }
}

template<class PermeabilityType>
void Foam::wells::generic<PermeabilityType>::sourceWell()
{
        // Access to current well cells and some convenient refs 
        const labelList& cells = wellSet_.toc();

        switch (operation_)
        {
            case well::PROD:
            {
                produceRate();
                break;
            }
            case well::INJE:
            {
                injectRate();
                break;
            }
        }

}

template<class PermeabilityType>
void Foam::wells::generic<PermeabilityType>::injectRate()
{
    const labelList& cells = wellSet_.toc();
    switch (mode_)
    {
        case well::LRATE:
        {
            // DISCUSSION: Pay close attention to "signs"
            // An injection rate goes to the rhs as a positive quantity
            // And production rates go there as negative quanitities
            // INFO: This strongly depends on the fact that a cell may belong
            // to "at most" one well
            forAll(cells, celli)
            {
                label cellID = cells[celli];
                // FIXME: Injection should be resisted by other fluids.. no??
                ewSource_[celli] = -targetRate().value()/mesh_.V()[cellID];
            }
            break;
        }
        case well::BHP:
        {
            forAll(cells, celli)
            {
                label cellID = cells[celli];
                ewSource_[celli] = (BHPData_.pi[celli]*kr1_[cellID]/mu1_[cellID])
                    *(BHPData_.BHP.value()-g_[2].value()*rho1_[cellID]*(BHPData_.DZ-mesh_.C()[cellID].z()));
                iwSource_[celli] = -BHPData_.pi[celli]*kr1_[cellID]/mu1_[cellID];
            }
            break;
        }
    }
}

template<class PermeabilityType>
void Foam::wells::generic<PermeabilityType>::produceRate()
{
    const labelList& cells = wellSet_.toc();
    switch (mode_)
    {
        case well::LRATE:
        {
            forAll(cells, celli)
            {
                label cellID = cells[celli];
                ewSource_[celli] = waterFrac_[cellID]*targetRate().value()
                                                    /mesh_.V()[cellID];
                eoSource_[celli] = (1.0-waterFrac_[cellID])*targetRate().value()
                                                    /mesh_.V()[cellID];
            }
            break;
        }
        case well::BHP:
        {
            forAll(cells, celli)
            {
                label cellID = cells[celli];
                Info << ewSource_ << endl;
                ewSource_[celli] = -(BHPData_.pi[celli]*kr1_[cellID]/mu1_[cellID])
                    *(BHPData_.BHP.value()-g_[2].value()*rho1_[cellID]*(BHPData_.DZ-mesh_.C()[cellID].z()));
                eoSource_[celli] = -(BHPData_.pi[celli]*kr2_[cellID]/mu2_[cellID])
                    *(BHPData_.BHP.value()-g_[2].value()*rho2_[cellID]*(BHPData_.DZ-mesh_.C()[cellID].z()));
                iwSource_[celli] = BHPData_.pi[celli]*kr1_[cellID]/mu1_[cellID];
                ioSource_[celli] = BHPData_.pi[celli]*kr2_[cellID]/mu2_[cellID];
            }
            break;
        }
    }
}

template<class PermeabilityType>
void Foam::wells::generic<PermeabilityType>::correct()
{

        ewSource_.setSize(wellSet_.toc().size());
        eoSource_.setSize(wellSet_.toc().size());
        iwSource_.setSize(wellSet_.toc().size());
        ioSource_.setSize(wellSet_.toc().size());
        switch (mode_)
        {
            case well::BHP:
            {
                calculateWellIndex();
                break;
            }
            case well::LRATE:
            {
                break;
            }
        }
        sourceWell();
}
