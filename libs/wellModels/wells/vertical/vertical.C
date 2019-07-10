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

#include "vertical.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class PermeabilityType>
Foam::wells::vertical<PermeabilityType>::vertical
(
    const word& name,
    const dictionary& wellsProperties,
    const fvMesh& mesh,
    const HCmixture<blackoilPhase>& mixture
)
:
    well(name, wellsProperties, mesh, mixture),
    p_(mesh.lookupObject<volScalarField>("p")),
    kr1_(mesh.lookupObject<volScalarField>(mixture.phase("wetting").name()+".kr")),
    kr2_(mesh.lookupObject<volScalarField>(mixture.phase("nonWetting").name()+".kr")),
    mu1_(mesh.lookupObject<volScalarField>(mixture.phase("wetting").name()+".mu")),
    mu2_(mesh.lookupObject<volScalarField>(mixture.phase("nonWetting").name()+".mu")),
    rho1_(mesh.lookupObject<volScalarField>(mixture.phase("wetting").name()+".rho")),
    rho2_(mesh.lookupObject<volScalarField>(mixture.phase("nonWetting").name()+".rho")),
    rFVF1_(mesh.lookupObject<volScalarField>(mixture.phase("wetting").name()+".rFVF")),
    rFVF2_(mesh.lookupObject<volScalarField>(mixture.phase("nonWetting").name()+".rFVF")),
    waterFrac_(mesh.lookupObject<volScalarField>("Fw")),
    g_
    (
        mesh.lookupObject<uniformDimensionedVectorField>("g")
    ),
    BHPData_(),
    weightCellRates_(wellDict_.lookupOrDefault<bool>("weightCellRates",false))
{
    // Conditionally initialize well control variables
    if (mode_ == well::LRATE)
    {
        flowRateData_ = flowRateModeVars(wellDict_);
        Info<< endl << "Well " << name_
            << " is operating in " << rateHandlingToWord(mode_) << " mode."
            << endl;
    } else if (mode_ == well::BHP)
    {
        BHPData_ = BHPModeVars(name_,wellDict_,mesh_);
        Info<< endl << "Well " << name_
            << " is operating in " << rateHandlingToWord(mode_) << " mode."
            << endl;
    } else {
        FatalErrorIn("Foam::well::well(const word& name, const dictionary& wellDict, const fvMesh& mesh)")
            << "Well control mode was not recognized; available conrol modes: 2(flowRate BHP);" 
            << exit(FatalError);
    }


    ewSource_.setSize(wellSet_.toc().size());
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

template<class PermeabilityType>
void Foam::wells::vertical<PermeabilityType>::calculateWellIndex()
{
    BHPData_.h.setSize(wellSet_.size());
    BHPData_.pi.setSize(wellSet_.size());
    // Const ref to absolute permeability field
    const PermeabilityType& k = mesh_.objectRegistry::lookupObject<PermeabilityType>("K");

    // First estimate cell sizes
    forAll(wellSet_.toc(), celli)
    {
        label cellID = wellSet_.toc()[celli];
        const labelList& cellEdges = mesh_.cellEdges()[cellID];
        scalarList dEdgesAlongZ(cellEdges.size());
        scalarList dEdgesAlongX(cellEdges.size());
        scalarList dEdgesAlongY(cellEdges.size());
        forAll(cellEdges, edgei){
            dEdgesAlongX[edgei] = mag(mesh_.edges()[cellEdges[edgei]].vec(mesh_.points()).x());
            dEdgesAlongY[edgei] = mag(mesh_.edges()[cellEdges[edgei]].vec(mesh_.points()).y());
            dEdgesAlongZ[edgei] = mag(mesh_.edges()[cellEdges[edgei]].vec(mesh_.points()).z());
        }
        BHPData_.h[celli][0] = gMax(dEdgesAlongX);
        BHPData_.h[celli][1] = gMax(dEdgesAlongY);
        BHPData_.h[celli][2] = gMax(dEdgesAlongZ);
    }


    // Estimate re
    if(wellDict_.found("re"))
    {
        // Get Equivalent radius for each cell from Dictionary
        wellDict_.lookup("re") >> BHPData_.re;
    } else
    {
        // If Equivalent radius is not explicitly provided, provide it
        estimateEquivRadius(k);
    }

    // Calculate well PI
    forAll(wellSet_.toc(), celli)
    {
        label cellID = wellSet_.toc()[celli];
        // Calculate sqrt(k11*k22)
        scalar Kfact = sqrtKK(k, cellID);

        // Then estimate cell PI
        BHPData_.pi[celli] = 2*mathematicalConstant::pi*Kfact*BHPData_.h[celli][2]
            /(log(BHPData_.re[celli]/radius_.value()) + BHPData_.skin);
    }

    // Debug code
    if(debug)
    {
        Info<< "Well " << name_ << " :" << nl
            << "H:" << tab << BHPData_.h << nl
            << "Re:" << tab << BHPData_.re << nl
            << "PI:" << tab << BHPData_.pi << endl;
    }
}

template<class PermeabilityType>
void Foam::wells::vertical<PermeabilityType>::sourceWell()
{
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
void Foam::wells::vertical<PermeabilityType>::injectRate()
{
    const labelList& cells = wellSet_.toc();
    switch (mode_)
    {
        case well::LRATE:
        {
            forAll(cells, celli)
            {
                label cellID = cells[celli];
                // FIXME: Forcing uniform injection into well
                ewSource_[celli] = -targetRate().value()*rFVF1_[cellID]/tV_; 
                ///mesh_.V()[cellID];
            }
            break;
        }
        case well::BHP:
        {
            if (switchedToBHP_) 
            {
                calculateWJMatrixCoeffs();
            } else {
                // Ordinary bhp-constant injection
                forAll(cells, celli)
                {
                    label cellID = cells[celli];
                    scalar WJ1 = BHPData_.pi[celli]*kr1_[cellID]*rFVF1_[cellID]/mu1_[cellID];
                    ewSource_[celli] = -WJ1*(BHPData_.BHP[0]-g_[2].value()*rho1_[cellID]*(BHPData_.DZ-mesh_.C()[cellID].z()));
                    iwSource_[celli] =  WJ1;
                }
            }
            break;
        }
    }
}

template<class PermeabilityType>
void Foam::wells::vertical<PermeabilityType>::produceRate()
{
    const labelList& cells = wellSet_.toc();
    switch (mode_)
    {
        case well::LRATE:
        {
            forAll(cells, celli)
            {
                label cellID = cells[celli];
                ewSource_[celli] = waterFrac_[cellID]*targetRate().value()*rFVF1_[cellID]/tV_;
                                                    ///mesh_.V()[cellID];
                eoSource_[celli] = (1.0-waterFrac_[cellID])*targetRate().value()*rFVF2_[cellID]/tV_;
                                                    ///mesh_.V()[cellID];
            }
            break;
        }
        case well::BHP:
        {
            forAll(cells, celli)
            {
                label cellID = cells[celli];
                scalar WJ1 = BHPData_.pi[celli]*kr1_[cellID]*rFVF1_[cellID]/mu1_[cellID]/mesh_.V()[cellID];
                scalar WJ2 = BHPData_.pi[celli]*kr2_[cellID]*rFVF2_[cellID]/mu2_[cellID]/mesh_.V()[cellID];
                ewSource_[celli] = WJ1*(BHPData_.BHP[0]-g_[2].value()*rho1_[cellID]*(BHPData_.DZ-mesh_.C()[cellID].z()));
                eoSource_[celli] = WJ2*(BHPData_.BHP[0]-g_[2].value()*rho2_[cellID]*(BHPData_.DZ-mesh_.C()[cellID].z()));
                iwSource_[celli] = -WJ1;
                ioSource_[celli] = -WJ2;
            }
            break;
        }
    }
}

template<class PermeabilityType>
void Foam::wells::vertical<PermeabilityType>::correct()
{

        ewSource_.setSize(wellSet_.toc().size(),0);
        eoSource_.setSize(wellSet_.toc().size(),0);
        iwSource_.setSize(wellSet_.toc().size(),0);
        ioSource_.setSize(wellSet_.toc().size(),0);
        switch (mode_)
        {
            case well::BHP:
            {
                // Handle Iso and Aniso permeabilities
                calculateWellIndex();
                break;
            }
            case well::LRATE:
            {
                // Calculate Dt until next rate change
                if(flowRateData_.flowRateSeries.empty()) break;
                forAll(flowRateData_.flowRateSeries, item)
                {
                    if((mesh_.time().value() != mesh_.time().startTime().value()) 
                        and 
                        (flowRateData_.flowRateSeries[item].first() >= mesh_.time().value()))
                    {
                        timeForDt_ = flowRateData_.flowRateSeries[item].first() - mesh_.time().value();
                        break;
                    }
                    break;
                }
                if (mesh_.time() > mesh_.time().startTime()) switchWellMode();
                break;
            }
        }
        sourceWell();
        BHPData_.setKeys();
}

template<class PermeabilityType>
void Foam::wells::vertical<PermeabilityType>::switchWellMode()
{
    // Do not switch if not allowed
    if(not(allowSwitchingModes_)) return;

    switch (mode_)
    {
        case well::BHP:
        {
            // Do nothing, this is the preferred way for now 
            break;
        }
        case well::LRATE:
        {
            // Switch to BHP mode but keep calculating new BHPs
            // according to flowRate at each timeStep
            BHPData_ = BHPModeVars(name_,wellDict_,mesh_,"switch");
            // Setup well kr model
            calculateWellIndex();
            switchedToBHP_ = true;
            //BHPData_.BHP.setSize(wellSet_.size());
            //updateBHP(0);
            mode_ = well::BHP;
            Info<< nl << "Well " << name_
                << " is switched to " << rateHandlingToWord(mode_) << " mode."
                << nl << endl;
            break;
        }
    }
}
template<class PermeabilityType>
void Foam::wells::vertical<PermeabilityType>::updateBHP()
{
    const labelList& cells = wellSet_.toc();
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    switch (operation_)
    {
        case well::INJE:
        {
            // Calculate BHP at first cell
            // Does this mean we end up with two different BHPs in parallel runs??
            scalar WJ1 = 0, WJ1p = 0, WJ1gammaDz = 0;
            forAll(cells, celli)
            {
                const label& cellID = cells[celli];
                WJ1 += BHPData_.pi[celli]*kr1_[cellID]*rFVF1_[cellID]/mu1_[cellID];
                WJ1p += p[cellID]*BHPData_.pi[celli]*kr1_[cellID]*rFVF1_[cellID]/mu1_[cellID];
                WJ1gammaDz += BHPData_.pi[celli]*kr1_[cellID]*rFVF1_[cellID]
                    *g_[2].value()*rho1_[cellID]*(BHPData_.DZ-mesh_.C()[cellID].z())/mu1_[cellID];
            }

            BHPData_.BHP[0] = (targetRate().value() + WJ1p - WJ1gammaDz)/WJ1;
            break;
        }
        case well::PROD:
        {
            // Do nothing
            break;
        }
    }
}
template<class PermeabilityType>
void Foam::wells::vertical<PermeabilityType>::calculateWJMatrixCoeffs()
{
    // Usefull references
    const volScalarField& p = mesh_.objectRegistry::lookupObject<volScalarField>("p");
    const labelList& cells = wellSet_.toc();

    /// Qi = WJi (BHP -pi -gammai.Dzi)

    // Calculate coeffs
    calculateWellIndex();
    scalarList WJ1(cells.size(), 0.0);
    scalarList WJ1p(cells.size(), 0.0);
    scalarList WJ1gammaDz(cells.size(), 0.0);
    scalarList WJ1coeffsForP(cells.size(), 0.0);
    forAll(cells, celli)
    {
        const label& cellID = cells[celli];
        WJ1[celli] = BHPData_.pi[celli]*kr1_[cellID]*rFVF1_[cellID]/mu1_[cellID];
        WJ1p[celli] = p[cellID]*WJ1[celli];
        WJ1gammaDz[celli] = BHPData_.pi[celli]*kr1_[cellID]*rFVF1_[cellID]
            *g_[2].value()*rho1_[cellID]*(BHPData_.DZ-mesh_.C()[cellID].z())/mu1_[cellID];
    }
    forAll(cells, celli)
    {
        const label& cellID = cells[celli];
        WJ1coeffsForP[celli] = WJ1[celli]*sum(WJ1p)/sum(WJ1)/mesh_.V()[cellID];
        WJ1coeffsForP[celli] -= WJ1[celli]*p[cellID]*WJ1[celli]/sum(WJ1)/mesh_.V()[cellID];
    }

    // The magic
    forAll(cells, celli)
    {
        const label& cellID = cells[celli];
        if(BHPData_.implicitWJp)
        {
            // Extremely dangerous because it meddles with top-level code for
            // matrices, you have been warned
            // Face addressing is painful to work with 
            //
            // Let's break it down
            // Assume I want to construct WJ.p
            // Assume i is owner of a shared face (index k) with a neighbor j
            // The coefficient multiplying p_j in row i of WJ 
            // is the (k)th member in upper() array
            // The coefficient multiplying p_i in row j of WJ 
            // is the (k)th member in lower() array
            
            // So, In the hope of this being stable,
            // and because it's the only option,
            // only pressure in cell neighbors in treated implicitly

            scalarField& lower = WJw_.lower();
            scalarField& upper = WJw_.upper();

            // Faces owned and neighbored by current cell
            labelList owned = findIndices(mesh_.owner(), cellID);
            labelList neighbored = findIndices(mesh_.neighbour(), cellID);

            if(debug > 3)
            Info<< "Cell " << cellID << " ownes " << owned << " and "
                << "neighbors " << neighbored <<  " faces." << endl;

            // loop through owned faces
            forAll(owned, owi)
            {

                // face neighbor
                const label& nei = mesh_.neighbour()[owned[owi]];
                label cellj = findIndex(cells, nei);
                if(debug > 3)
                Info<< "The neighbor " << nei << " is at position "
                    << cellj << " in well cells." << endl;
                if ( cellj >= 0)
                {
                    // If the neighboring cell is also a well cell
                    lower[owned[owi]] = WJ1[cellj]*WJ1[celli]/sum(WJ1)/mesh_.V()[cellID];
                    Info << "Subtracting from explicit coeffs at cell " << cellj << " in  well."
                    << endl;
                    WJ1coeffsForP[cellj] -= WJ1[cellj]*p[cells[cellj]]*WJ1[cellj]/sum(WJ1)/mesh_.V()[cells[cellj]];
                } else {
                    lower[owned[owi]] = 0.0;
                }
            }

            
            // loop through neighbored faces
            forAll(neighbored, nei)
            {

                // face neighbor
                const label& owi = mesh_.owner()[neighbored[nei]];
                label cellj = findIndex(cells, owi);
                if(debug > 3)
                Info<< "The owner " << owi << " is at position "
                    << cellj << " in well cells." << endl;
                if ( cellj >= 0)
                {
                    // If the neighboring cell is also a well cell
                    upper[neighbored[nei]] = WJ1[cellj]*WJ1[celli]/sum(WJ1)/mesh_.V()[cellID];
                    Info << "Subtracting from explicit coeffs at cell " << cellj << " in  well."
                    << endl;
                    WJ1coeffsForP[cellj] -= WJ1[cellj]*p[cells[cellj]]*WJ1[cellj]/sum(WJ1)/mesh_.V()[cells[cellj]];
                } else {
                    upper[neighbored[nei]] = 0.0;
                }
            }
            if(debug > 3)
            Info << endl;

        ewSource_[celli] = -WJ1[celli]*(targetRate().value()+sum(WJ1gammaDz))/sum(WJ1)/mesh_.V()[cellID]
            + WJ1[celli]*g_[2].value()*rho1_[cellID]*(BHPData_.DZ-mesh_.C()[cellID].z())/mesh_.V()[cellID]
            - WJ1coeffsForP[celli];
        iwSource_[celli] = WJ1[celli]/mesh_.V()[cellID] 
            + WJ1[celli]*p[cellID]*WJ1[celli]/sum(WJ1)/mesh_.V()[cellID];

        } else {
            // Only adds to diagonal
            // Sparseness insured but not stability
            forAll(cells, cellii)
            {
                if (cellii == celli) continue;
                const label& j = cells[cellii];
                ewSource_[celli] += -WJ1[cellii]*p[j]/sum(WJ1)/mesh_.V()[celli];
            }
            iwSource_[celli] += -WJ1[celli]*p[cellID]/sum(WJ1)/mesh_.V()[celli];
        }

    }

    if(debug < 4) return;

    Info<< nl << "##################" << nl;
    List<List<scalar> > mat(WJw_.diag().size(), List<scalar>(WJw_.diag().size(), 0.0));
    
    for(label k=0; k<WJw_.diag().size(); k++) mat[k][k] = WJw_.diag()[k];
    for(label k=0; k<WJw_.lduAddr().lowerAddr().size(); k++)
    {
        label i = WJw_.lduAddr().lowerAddr()[k];
        label j = WJw_.lduAddr().upperAddr()[k];
        mat[i][j] = WJw_.upper()[k];
        mat[j][i] = WJw_.lower()[k];
    }
    
    for(label row=0; row<mat.size(); row++)
    {
        forAll(mat[row], item)
        {
            Info << mat[row][item] << ",";
        }
        Info << nl;
    }
}


template<class PermeabilityType>
bool Foam::wells::vertical<PermeabilityType>::writeData(Ostream& os) const
{
    os.writeKeyword("orientation");
    os <<  type_ << ";" << endl;

    os.writeKeyword("rateMode");
    os <<  rateHandlingToWord(mode_) << ";" << endl;

    os.writeKeyword("operationMode");
    os <<  operationHandlingToWord(operation_) << ";" << endl;

    os.writeKeyword("radius");
    os <<  radius_ << ";" << endl << nl;

    os.writeKeyword("cells");
    os <<  cellIDs() << ";" << endl << nl;

    os.writeKeyword("totalVolume");
    os <<  tV_ << ";" << endl;

    os.writeKeyword("targetRate");
    if (mode_ == well::LRATE or switchedToBHP_)
    {
        os <<  targetRate() << ";" << endl;
    } else {
        os <<  0.0 << ";" << endl;
    }

    os.writeKeyword("deltaTForRateChange");
    os <<  timeForDt_ << ";" << endl;

    os.writeKeyword("allowSwitchingModes");
    os <<  allowSwitchingModes_ << ";" << endl;

    os.writeKeyword("switchedToBHP");
    os <<  switchedToBHP_ << ";" << endl;

    os.writeKeyword("BHPData");
    //BHPData_.setKeys();
    os <<  BHPData_ << endl;

    return os.good();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //
