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

Namespace
    Foam::capillaryPressureModels

Description
    A namespace for capillarity models.

Class
    Foam::capillaryPressureModelBase

Description
    An abstract class template for capillarity models.

SourceFiles
    capillaryPressureModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef capillaryPressureModelBase_H
#define capillaryPressureModelBase_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"
#include "volFieldsFwd.H"
#include "krPcCompileCalculations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<int nPhases>
using cmptT = typename pcTableReturn<nPhases>::cmptType;

/*---------------------------------------------------------------------------*\
                         Class capillaryPressureModelBase Declaration
\*---------------------------------------------------------------------------*/

template<int nPhases>
class capillaryPressureModelBase
{

protected:

    // Protected data

        //- Model Name
        word name_;

        //- Const-ref to Hosting dict
        const dictionary& transportProperties_;

        //- Const-ref to mesh
        const fvMesh& mesh_;

        //- Model capillary Pressure
        typename pcTableReturn<nPhases>::type pc_;

        //- Model capillary Pressure derivatives
        typename pcTableReturn<nPhases>::type dpcdS_;

        //- Disallow copy construct
        capillaryPressureModelBase(const capillaryPressureModelBase&);

        //- Disallow default bitwise assignment
        void operator=(const capillaryPressureModelBase&);

public:

    // Class typeName
    TypeName("capillaryPressureModelBase");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable ( 	  	
        autoPtr,
		capillaryPressureModelBase,
		dictionary,
		(
        const word& name,
        const dictionary& transportProperties,
        const fvMesh& mesh
        ),
        (name, transportProperties, mesh)
	); 		

    // Constructors

        //- Construct from components
        capillaryPressureModelBase
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

    // Selectors

        //- Return a reference to the selected Pc model
        static autoPtr<capillaryPressureModelBase<nPhases> > New
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~capillaryPressureModelBase() {};

    // Public Member Functions

        //- Return model hosting dictionary
        const dictionary& transportProperties() const
        {
            return transportProperties_;
        }

        //- Return capillarity component
        tmp<volScalarField> pc(int phaseIndex = 0)
        {
            return pcTableReturn<nPhases>::pcComponent(pc_, phaseIndex);
        }

        //- Return capillarity derivative
        tmp<volScalarField> dpcdS(int phaseIndex = 0)
        {
            return pcTableReturn<nPhases>::pcComponent(dpcdS_, phaseIndex);
        }

        //- Return Pc table as a hole (For debugging) by ref
        const typename pcTableReturn<nPhases>::type& pcTable() const
        {
            return pc_;
        }

        //- Return dPc table as a hole (For debugging) by ref
        const typename pcTableReturn<nPhases>::type& dpcdSTable() const
        {
            return dpcdS_;
        }

        //- Correct capillarity fields
        virtual void correct() = 0;

        //- Return associated phase names as a list
        virtual const wordList phaseNames() const = 0;

        //- Return canonical phase name
        virtual const word canonicalPhase() const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Macro for RTS-enabled class initiated from this class template
#define makeCapillaryPressureModelBase(nPhases)                             \
                                                                               \
defineNamedTemplateTypeNameAndDebug(base##nPhases##PhasesPcModel, 0);          \
defineTemplateRunTimeSelectionTable(base##nPhases##PhasesPcModel, dictionary); \

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "capillaryPressureModelBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
