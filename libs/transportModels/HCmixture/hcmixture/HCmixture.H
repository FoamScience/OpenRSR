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

Class
    Foam::HCmixture

Description
    A Generic templated fluid HCmixture class to handle mixture-related 
    properties of fluids.

SourceFiles
    HCmixture.C

\*---------------------------------------------------------------------------*/

#ifndef HCmixture_H
#define HCmixture_H

#include "bubblePointModel.H"
#include "regIOobject.H"
#include "volFields.H"
#include "List.H"
#include "phase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                         Class fluidPhase Declaration
\*---------------------------------------------------------------------------*/

template<class FluidModel>
class HCmixture
    : public regIOobject
{

protected:

    // Protected Data Members
        
        //- Mixture name
        word name_;

        // Phase Dict
        dictionary HCmixtureDict_;

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Phase names
        wordList phaseNames_;

        //- Wettability Status
        //  Options: wetting, nonWetting, intermediateWetting
        wordList wettingStatus_;

        //- List of refs to phase objects
        List<autoPtr<FluidModel> > phaseList_;

        //- Pointer to bubble point model
        autoPtr<bubblePointModel> pbModel_;
public:

    // Constructors

        //- Construct from components
        HCmixture
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

    // Selectors

        //- Return a reference to HCmixture
        static autoPtr<HCmixture<FluidModel> > New
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~HCmixture(){}

    // Public Member Functions
    
        //- Return HCmixture name
        const word& name() const
        {
            return name_;
        }
        
        //- Return const ref to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return const ref to phase names
        const wordList& phaseNames() const
        {
            return phaseNames_;
        }

        //- Return phases list
        const List<autoPtr<FluidModel> > phases() const
        {
            return phaseList_;
        }

        //- Return ref to phase by wettingStatus
        FluidModel& phase(word status)
        {
            forAll(wettingStatus_,phasei){
                if (wettingStatus_[phasei] == status)
                    return phaseList_[phasei]();
            }
            FatalErrorIn("Foam::FluidModel& Foam::HCmixture::phase(word status)")
                << "No " << status << " phase set." << nl << exit(FatalError);
        }

        //- Return ref to phase by wettingStatus
        const FluidModel& phase(word status) const
        {
            forAll(wettingStatus_,phasei){
                if (wettingStatus_[phasei] == status)
                    return phaseList_[phasei]();
            }
            FatalErrorIn("Foam::FluidModel& Foam::HCmixture::phase(word status)")
                << "No " << status << " phase set." << nl << exit(FatalError);
        }

        void correct()
        {
            // Correct phases
            forAll(phaseList_, phasei)
            {
                phaseList_[phasei]->correct();
            }

            // Correct Pb model
            pbModel_->correct();
        }

        bool writeData(Ostream&) const 
        {
            // Not implemented yet
            return 1;
        }

    // Public Member Operators    

        //- indexing operator [] to retrieve phases
        inline FluidModel& operator[](label i);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HCmixture.C"

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
