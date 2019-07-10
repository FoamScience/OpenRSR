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

\*---------------------------------------------------------------------------*/

#include "wellFlowRate.H"
#include "volFields.H"
#include "dictionary.H"
#include "foamTime.H"
#include "HCmixture.H"
#include "blackoilPhase.H"
#include "well.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wellFlowRate, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wellFlowRate::wellFlowRate
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    meshName_(fvMesh::defaultRegion),
    active_(true),
    log_(false),
    wellSet_(),
    wellFlowRateFilePtr_(NULL)
{
    // Check if the available mesh is an objectRegistry otherise deactivate
    if (!isA<objectRegistry>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "wellFlowRate::wellFlowRate"
            "(const objectRegistry& obr, const dictionary& dict)"
        )   << "No objectRegistry available, deactivating."
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wellFlowRate::~wellFlowRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wellFlowRate::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", false);
        dict.lookup("wells") >> wellSet_;
        meshName_ = dict.lookupOrDefault<word>("region", fvMesh::defaultRegion);
    }
}


void Foam::wellFlowRate::makeFile()
{
    // Create the wellFlowRate file if not already created
    if (wellFlowRateFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating wellFlowRate file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName wellFlowRateDir;
            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                wellFlowRateDir =
                    obr_.time().path()/".."/name_/"0";
            }
            else
            {
                wellFlowRateDir =
                    obr_.time().path()/name_/"0";
            }

            // Create directory if does not exist.
            mkDir(wellFlowRateDir);

            // Open new file at start up
            wellFlowRateFilePtr_.reset
            (
                new OFstream(wellFlowRateDir/(type() + ".dat"))
            );

            // Add headers to output data
            writeFileHeader();
        }
    }
}


void Foam::wellFlowRate::writeFileHeader()
{
    if (wellFlowRateFilePtr_.valid())
    {
        wellFlowRateFilePtr_()
            << "# Time" << tab << "well" << tab << "waterProdRate" 
            << tab << "waterInjRate" << tab << "oilProdRate"
            << endl;
    }
}


void Foam::wellFlowRate::execute()
{
    // Do nothing - only valid on write
}


void Foam::wellFlowRate::end()
{
    // Do nothing - only valid on write
}


void Foam::wellFlowRate::write()
{
    if (active_)
    {
        // Create the wellFlowRate file if not already created
        makeFile();


        // We can recalculate everything here on per well basis
        // but for efficiency & to reveal bugs in well models:
        // Get references to Well Sources coefficients (Explicit-Implicit)
        // in solver equations from well management model
        // // // // // 
        // The following is strongly based on the fact that every cell in 
        // the mesh belongs to --at most-- ONE well.

        // Access the pressure field
        const volScalarField& p = obr_.lookupObject<volScalarField>("p");
        const fvMesh& mesh = p.mesh();

        // So, access to production/injection coeffs
        const volScalarField& eoSource = mesh.lookupObject<volScalarField>("oilExpSource");
        const volScalarField& ewSource = mesh.lookupObject<volScalarField>("waterExpSource");
        const volScalarField& ioSource = mesh.lookupObject<volScalarField>("oilImpSource");
        const volScalarField& iwSource = mesh.lookupObject<volScalarField>("waterImpSource");

        volScalarField wRate(ewSource+iwSource*p);
        volScalarField oRate(eoSource+ioSource*p);

        // Now treat each well depending on its mode of operation :)
        forAll(wellSet_, welli)
        {
            word wellName = wellSet_[welli];

            scalar oR(0.0);
            scalar wR(0.0);

            if (mesh.foundObject<well>(wellName))
            {
                const well& thisWell = mesh.lookupObject<well>(wellName);
                const labelList cells = thisWell.cellIDs();
                const scalar tV = thisWell.totalVolume();

                forAll(cells, celli)
                {
                    const label& cellID = cells[celli];
                    wR += wRate[cellID]*mesh.V()[cellID];
                    oR += oRate[cellID]*mesh.V()[cellID];
                }

                if (log_)
                {
                    Pout<< "wellFlowRate: " 
                        << wellName << tab << wR
                        << tab << oR << endl;
                }
                // Write to file
                if (Pstream::master())
                {
                    wellFlowRateFilePtr_() << mesh.time().value() << tab
                        << wellName << tab << wR 
                        << tab << oR << endl;

                }
            } else {
                Info << name_ <<  ": Well " << wellName << " not found" << endl;
            }
        }
    }
}




// ************************************************************************* //
