/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::wellFlowRate

Description
    Calculates production & injection rates for selected wells

    Data written to the file \<startTimeDir\>/wellFlowRate.dat


SourceFiles
    wellFlowRate.C
    IOwellFlowRate.H

\*---------------------------------------------------------------------------*/

#ifndef wellFlowRate_H
#define wellFlowRate_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "OFstream.H"
#include "Switch.H"
#include "pointFieldFwd.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class fvMesh;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class wellFlowRate Declaration
\*---------------------------------------------------------------------------*/

class wellFlowRate
{
protected:

    // Protected data

        //- Name of this set of well flowrates
        //  Also used as the name of the output directory.
        word name_;

        //- Object registery, the mesh
        const objectRegistry& obr_;

        //- Mesh region name
        word meshName_;

        //- on/off switch
        bool active_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        //- Wells in question
        wordList wellSet_;

        //-Ptr to output file
        autoPtr<OFstream> wellFlowRateFilePtr_;


    // Private Member Functions

        //- If the output file has not been created create it
        void makeFile();

        //- Disallow default bitwise copy construct
        wellFlowRate(const wellFlowRate&);

        //- Disallow default bitwise assignment
        void operator=(const wellFlowRate&);

        //- Output file header information
        virtual void writeFileHeader();


public:

    //- Runtime type information
    TypeName("wellFlowRate");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        wellFlowRate
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

    //- Destructor
    virtual ~wellFlowRate();


    // Member Functions

        //- Return name of the set of wells flowrates
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the well data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Write the wellFlowRate
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "wellFlowRateTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
