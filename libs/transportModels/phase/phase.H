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
    Foam::phase

Description
    A Generic fluid phase class which registers created phases to the 
    mesh's objectRegistry.

SourceFiles
    phase.C

\*---------------------------------------------------------------------------*/

#ifndef phase_H
#define phase_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                         Class fluidPhase Declaration
\*---------------------------------------------------------------------------*/

class phase
{

protected:

    // Protected Data Members
        
        //- Phase name
        word name_;

        // Phase Dict
        dictionary phaseDict_;

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Phase velocity
        volVectorField U_;

        //- Phase saturation in a mixture
        volScalarField alpha_;

public:

    // Constructors

        //- construct from components
        phase
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

    // Selectors

        ////- Return a reference to phase
        //static autoPtr<phase> New
        //(
        //    const word& name,
        //    const dictionary& transportProperties,
        //    const fvMesh& mesh
        //);

    //- Destructor
    virtual ~phase(){}

    // Public Member Functions
    
        //- Return phase name
        const word& name() const
        {
            return name_;
        }
        
        //- Return const ref to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return const ref to velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return ref to velocity
        volVectorField& U()
        {
            return U_;
        }

        //- Return const ref to saturation
        const volScalarField& alpha() const
        {
            return alpha_;
        }

        //- Return ref to saturation
        volScalarField& alpha()
        {
            return alpha_;
        }

        // Correct phase properties
        virtual void correct() = 0;

        // Write phase Data
        bool writeData(Ostream& os) const 
        {
            return true;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
