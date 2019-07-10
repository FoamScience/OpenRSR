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
    Foam::relativePermeabilityModels::krBrooksCorey

Description
    Standard two-phases Brooks-Corey relative permeability model:
    \f[
        k_{rw} = k_{rw, max}S_{w, norm}^n
    \f]
    \f[
        k_{rn} = k_{rn, max}(1-S_{w, norm})^n
    \f]

    Where 
    \f[
        S_{w, norm} = \frac{S_w, S_{w, irr}}{1-S_{w, irr}-S_{n, irr}}
    \f]
    \f[
        S_{w, norm} = 0\ if\ S_w<S_{w, irr}
    \f]
    \f[
        S_{w, norm} = 1\ if\ S_w>1-S_{n, irr}
    \f]

    \vartable
        k_r             | Relative permeability
        w,n             | Wetting and NonWetting indices respectively
        S_{w, norm}     | Normalized wetting-phase saturation
        S_{w, irr}      | Irreducible wetting phase saturation (`Swmin` in `transportProperties`)
        S_{n, irr}      | Irreducible wetting phase saturation (`Swmax = 1-Sn_irr` in `transportProperties`)
    \endtable

SourceFiles
    krBrooksCorey.C

\*---------------------------------------------------------------------------*/

#ifndef krBrooksCorey_H
#define krBrooksCorey_H

#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------*\
        Class krBrooksCorey Declaration
\*---------------------------------------------------------------------------*/

class krBrooksCorey
:
    public relativePermeabilityModel
{
    //- Private Data Members

        //- Coefficients dictionary
        dictionary coeffsDict_;
        
        //- Residual wetting saturation
        volScalarField Swc_;
        
        //- 1-Residual nonwetting saturation
        volScalarField Sor_;

        //- Model exponent
        volScalarField n_;

        //- Normalized wetting saturation 
        volScalarField Se_;

        //- Max permeabilities
        volScalarField kr1Max_,kr2Max_;
        
public:

    //- Runtime type information
    TypeName("BrooksCorey");

    // Constructors

    //- Construct from components
    krBrooksCorey
    (
        const word& name,
        const dictionary& transportProperties,
        const phase& phase1,
        const phase& phase2
    );

    //- Destructor
    ~krBrooksCorey()
        {}

    // Member Functions

    //- Correct the relative permeabilities
    void correct();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
