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
    Foam::wells::BHPModeVars

Description
    A structure to hold BHP-related well data

\*---------------------------------------------------------------------------*/

#ifndef BHPModeVars_H
#define BHPModeVars_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wells
{

struct BHPModeVars
:
    public dictionary
{
    //- Well PI
    scalarList pi;

    //- Cells approx. sizes
    List<vector> h;

    //- Datum depth
    scalar DZ;

    //- BHP datum depth
    //- BHP at each cell if  injector converted to BHP-mode
    scalarList BHP;

    //- well equivalent radius
    scalarList re;

    //- well Skin
    scalar skin;

    //- Switch to update BHP at each timeStep
    bool keepUpdatingBHP;

    //- Turn on/off fully implicit treatment for pressure in all 
    //  well cells while injecting with BHP-operated wells
    bool implicitWJp;

    // Constructors

        //- Default Constructor
        BHPModeVars()
        :
        dictionary(),
        pi(1),
        h(1),
        DZ(0.0),
        BHP(1),
        re(1),
        skin(0.0),
        keepUpdatingBHP(false),
        implicitWJp(false)
        {
            createKeys();
        }

        //- Construct from Components
        BHPModeVars
        (   
            const word& name,
            const dictionary& wellDict,
            const fvMesh& mesh
        )
        :
        dictionary(),
        pi(1),
        h(1),
        DZ(readScalar(wellDict.lookup("datumDepth"))),
        BHP(1),
        re(1),
        skin(readScalar(wellDict.lookup("skin"))),
        keepUpdatingBHP(false),
        implicitWJp(false)
        {
            wellDict.lookup("BHP") >> BHP[0];
            createKeys();
        }

        //- Constructor for use ONLY when switching modes
        BHPModeVars
        (   
            const word& name,
            const dictionary& wellDict,
            const fvMesh& mesh,
            const word&
        )
        :
        dictionary(),
        pi(1),
        h(1),
        DZ(readScalar(wellDict.lookup("datumDepth"))),
        BHP(1),
        re(1),
        skin(readScalar(wellDict.lookup("skin"))),
        keepUpdatingBHP(true),
        implicitWJp(wellDict.lookupOrDefault<bool>("implicitWJpInNeighbors",false))
        {
            createKeys();
        }

        //- Populate the BHP dictionary with enteries
        void createKeys()
        {
            add("PI", pi);
            add("h", h);
            add("DZ", DZ);
            add("BHP", BHP);
            add("re", re);
            add("skin", skin);
            add("keepUpdatingBHP", keepUpdatingBHP);
            add("implicitWJp", implicitWJp);
        }

        //- Populate the BHP dictionary with enteries
        void setKeys()
        {
            set("PI", pi);
            set("h", h);
            set("DZ", DZ);
            set("BHP", BHP);
            set("re", re);
            set("skin", skin);
            set("keepUpdatingBHP", keepUpdatingBHP);
            set("implicitWJp", implicitWJp);
        }
};

}

}
#endif
