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
    Foam::wells::vertical

Description
    There should be an obvious separation for vertical, horizontal
    and inclined wells in the way a well model is created.

    Support for horizontal well models is planned but not yet thought
    through; i.e. It's safe to make this class change or vanish at will as
    long as well base class loses its virtuality.

    DISCLAIMER: Class not tested in parallel env.

SourceFiles
    vertical.C

\*---------------------------------------------------------------------------*/

#ifndef vertical_H
#define vertical_H

#include "well.H"
#include "volFields.H"
#include "HCmixture.H"
#include "cellToFace.H"
#include "BHPModeVars.H"
#include "blackoilPhase.H"
#include "relativePermeabilityModel.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wells
{

/*---------------------------------------------------------------------------*\
                  Class vertical Declaration
\*---------------------------------------------------------------------------*/

template<class PermeabilityType>
class vertical
:
    public well
{

protected:

        //- Const ref to pressure field
        const volScalarField& p_;

        //- Refs to relative permeabilities
        const volScalarField& kr1_;
        const volScalarField& kr2_;

        //- Refs to phase properties
        const volScalarField& mu1_;
        const volScalarField& mu2_;
        const volScalarField& rho1_;
        const volScalarField& rho2_;
        const volScalarField& rFVF1_;
        const volScalarField& rFVF2_;

        //- Ref to water fractional flow;
        //  Should be calculated once in solver
        //  FIXME: We are interested in knowing Fw only in Well cells if LRATE is
        //  specified, aren't we?
        //  Why calculate it everywhere then??
        //  DISCUSSION: Complications of making this abstract
        //  Qw/(Qw+Qo) should use Qi = -(K kri A/mu)*grad(p)
        //  and if grad(p) is volScalarField we're good to go
        const volScalarField& waterFrac_;

        //- Gravity field
        uniformDimensionedVectorField g_;

        //- BHP data
        BHPModeVars BHPData_;

        //- Should we weight injection into cell with permeability?
        bool weightCellRates_;

        //- Private Member Functions

            //- Calculate sqrt(k11*k22) for PI
            scalar sqrtKK(const PermeabilityType& k, const label& cellID);

public:

    //- Runtime type information
    TypeName("vertical");

    // Constructors

    //- Construct from components
    vertical
    (
        const word& name,
        const dictionary& wellDict,
        const fvMesh& mesh,
        const HCmixture<blackoilPhase>& mixture
    );

    //- Destructor
    ~vertical() {}

    // Member Functions
    
        //- Return BHP
        const dimensionedScalar& bhp() const {
            return BHPData_.BHP;
        }

        //- Return Well Productivity Index
        const scalarList& pi() const {
            return BHPData_.pi;
        }

        //- Return Well datum depth
        const scalar& dz() const {
            return BHPData_.DZ;
        }

        //- Estimate well equivalent radius (Must be specialized)
        void estimateEquivRadius(const PermeabilityType& k);

        //- Return coefficients for injection cells (Msut be specialized)
        scalarList correctInjectionIntoCells(const PermeabilityType& k);

        //- Calculate well index
        void calculateWellIndex();

        //- Calculate well source
        void sourceWell();

        //- Correct coeffs for fixed-rate water injection
        void injectRate();

        //- Correct coeffs for fixed-rate production
        void produceRate();

        //- Correct coeffs for fixed-rate production
        void correct();

        //- Switch rate-imposed and BHP-imposed modes in a well 
        void switchWellMode();

        //- Calculate WJ.p coeffs for implicit p in wells
        void calculateWJMatrixCoeffs();

        //- Update BHP at first cell in well
        void updateBHP();

        //- Write well data
        bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wells

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
