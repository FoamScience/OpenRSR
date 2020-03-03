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
    Foam::twoPhaseWellModels::peaceman

Description
    A class implementing a "workaround" for removing singularity when
    solving for wells BHP; Taking the most simple case. The volumetric
    flow rate of a phase $i$ at the well is expressed as follows:

    \f[
        Q_{i} = J\frac{k_{ri}}{\mu_{i}}[BHP - p_i -\gamma_i(z_{BH}-z)]
    \f]
    \f[
        J = \frac{2\phi h K}{ln(\frac{r_e}{r_w})+s}
    \f]
    \f[
        r_e = 0.14\sqrt{h_1^2+h_2^2}
    \f]

    Where
    \vartable
        J               | Well Index
        r_e             | Well equivalent radius
        k_{ri}          | Relative permeability of phase i in a well cell
        \mu_i           | Phase dynamic viscosity in a well cell
        z_{BH}          | The depth of the Bottom-Hole point
        z               | Cell Depth
        \gamma_i        | Fluid weight
    \endtable

SourceFiles
    peaceman.C

\*---------------------------------------------------------------------------*/

#ifndef peaceman_H
#define peaceman_H

#include "wellBaseAndModel.H"
#include "relativePermeabilityModelBase.H"
#include "capillaryPressureModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
//namespace twoPhaseWellModels
//{

/*---------------------------------------------------------------------------*\
                  Class peaceman Declaration
\*---------------------------------------------------------------------------*/

template<class KType, class MuRhoType>
class peaceman
: public wellModelBase<KType, 2>
{

protected:

    // Protected Data Members 

        //- Ref to Kr Model
        const relativePermeabilityModelBase<2>& krModel_;

        //- Ref to pc Model
        const capillaryPressureModelBase<2>& pcModel_;

        //- Canonical Phase viscosity
        const MuRhoType& canonicalMu_;

        //- Non-Canonical Phase viscosity
        const MuRhoType& nonCanonicalMu_;

        //- Canonical Phase density
        const MuRhoType& canonicalRho_;

        //- Non-Canonical Phase density
        const MuRhoType& nonCanonicalRho_;

    // Protected Member Functions

        //- Calculate FractionalFlow In well cells
        scalarList calculateFractionalFlow(const wellBase<KType, 2>& well) const;

        //- Calculate the ratio totalRate/cellRate
        scalar calculateCellRateRatio
        (
            const label& cellID,
            const wellBase<Iso, 2>& well,
            bool forCanonical
        ) const;

        //- Estimate Well Cell sizes
        List<vector> estimateCellSizes(const wellBase<KType, 2>& well) const;

        //- Calculate Well Equivalent Radius
        scalarList estimateEquivRadius(const wellBase<KType, 2>& well) const;

        //- Calculate Well Productivity Index
        scalarList calculateWellPI(const wellBase<KType, 2>& well) const;

public:

    //- Runtime type information
    TypeName("peaceman");

    // Declare run-time constructor selection table
    //declareRunTimeSelectionTable
    //(
    //    autoPtr,
    //    wellModelBase,
    //    Dictionary,
    //    (
    //        const word& name,
    //        const dictionary& wellProperties,
    //        const fvMesh& mesh
    //    ),
    //    (name, wellProperties, mesh)
    //);

    // Constructors

        //- Construct from components
        peaceman
        (
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh
        );

    //- Destructor
    ~peaceman() {}

    // Member Functions

        //- Correct well model
        virtual void correct();

        //- Calculate Well source
        virtual void operator()(const word& wellName) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//}; // End namespace twoPhaseWellModels

// Macro for RTS-enabled class initiated from this class template
#define makePeacemanModel(KType, MuRhoType)                                     \
                                                                                \
defineNamedTemplateTypeNameAndDebug(derived##KType####MuRhoType##Peaceman, 0);    \
addToRunTimeSelectionTable                                                      \
(                                                                               \
    base2##KType##WellModel,                                              \
    derived##KType####MuRhoType##Peaceman,                                      \
    Dictionary                                                                  \
);


} // End namespace Foam

#ifdef NoRepository
#   include "peaceman.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
//makeWellModelDerived(Peaceman, KType, MuRhoType);                               \
//defineTemplateRunTimeSelectionTable                                             \
//(                                                                               \
//    derived##KType####MuRhoType##Peaceman,                                      \
//    Dictionary                                                                  \
//);                                                                              \
