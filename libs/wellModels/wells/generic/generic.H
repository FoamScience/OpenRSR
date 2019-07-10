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
    Foam::wells::generic

Description
    Generic arbitrary-direction, (to-be) multi-branch well modeling in 
    unstructured mesh.

SourceFiles
    generic.C

\*---------------------------------------------------------------------------*/

#ifndef generic_H
#define generic_H

#include "fvCFD.H"
#include "vertical.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wells
{

/*---------------------------------------------------------------------------*\
                  Class generic Declaration
\*---------------------------------------------------------------------------*/

template<class PermeabilityType>
class generic
:
    public well
{
    // Private Member Variables
        //- Refs to relative permeabilities
        const volScalarField& kr1_;
        const volScalarField& kr2_;

        //- Refs to phase properties
        const volScalarField& mu1_;
        const volScalarField& mu2_;
        const volScalarField& rho1_;
        const volScalarField& rho2_;

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

        //- Name of global transmissibility field
        word TfName_;

        //- Name of global mu field
        word muName_;

        //- Const-ref to global pressure field
        const volScalarField& p_;
 
        //- Const-ref to global transmissibility field
        const surfaceScalarField& Tf_;

        //- Const-ref to global viscosity field
        const volScalarField& mu_;

        //- Mesh subset
        fvMeshSubset meshSubset_;

        //- Cell map from origin to subset
        labelList cellmap_;

        //- Face map from origin to subset
        labelList facemap_;

        //- A subset of original mesh
        const fvMesh& submesh_;

        //- Test pressure field
        autoPtr<volScalarField> tP_;

        //- Test flowrate per volume unite field
        autoPtr<volScalarField> tQ_;

        //- Local transmissiblity
        autoPtr<surfaceScalarField> tTf_;

    // Private Member Function
        
public:

    //- Runtime type information
    TypeName("generic");

    // Constructors

        //- Construct from components
        generic
        (
            const word& name,
            const dictionary& wellDict,
            const fvMesh& mesh,
            const HCmixture<blackoilPhase>& mixture
        );

        //- Destructor
        ~generic() {}

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

        //- Write well data
        bool writeData(Ostream&) const
        {
            // Not implemented yet
            return 1;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wells

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
