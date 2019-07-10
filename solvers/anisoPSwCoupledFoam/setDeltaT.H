//// Timestepping Saturation equation
//volScalarField dSdt = mag(fvc::div(phiw)+ewSource+iwSource*p)/porosity;
//scalar dtForS = dSmax/(max(dSdt).value()+SMALL);
//
//scalar CFLUse = -1;
//scalar maxDeltaTFact = -1;
//#include "impesCoatsNo.H"
//
//scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
//
//// Set the timeStep
////runTime.setDeltaT
////(
////    min(
////        dtForS,
////        min(
////            deltaTFact*runTime.deltaTValue(),
////            maxDeltaT
////        )
////    )
////);
//
//// Homogeneous DeltaT in parallel runs
//// Would it be harmful if we had a different timeStep in each processor region??
//scalarList gatheredDeltaT(Pstream::nProcs());
//gatheredDeltaT[Pstream::myProcNo()] =
//(
//    min(
//        dtForS,
//        min(
//            deltaTFact*runTime.deltaTValue(),
//            maxDeltaT
//        )
//    )
//);
//Pstream::gatherList(gatheredDeltaT);
//
//runTime.setDeltaT(min(gatheredDeltaT));
//
//Info<< "deltaT = " <<  runTime.deltaTValue() << endl;

// ************************************************************************* //
// Timestepping Saturation equation
scalarField dSdt= mag(fvc::div(phiw)+ewSource+iwSource*p)().internalField()/porosity.internalField();
scalar dtForS = dSmax/(gMax(dSdt)+SMALL);

scalar CFLUse = -1;
scalar maxDeltaTFact = -1;
///// New impesCoatsNo for anisotropy
#include "impesCoatsNo.H"

scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
scalar cflFactor = runTime.controlDict().lookupOrDefault<scalar>("cflFactor", 2.0);
scalar minimalDeltaT = runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", 100);

// wells changing imposed rate
scalar wellDt = wModel->nextTimeStep();

// Set the timeStep
runTime.setDeltaT
(
        min(
            wellDt,
            min(
                dtForS,
                cflFactor*min(
                    deltaTFact*runTime.deltaTValue(),
                    maxDeltaT
                )
            )
        )
);

if (runTime.deltaTValue() < minimalDeltaT)
{
    Info<< endl << "Failed because of small deltaT." << endl;
    runTime.writeAndEnd();
}

Info<< "deltaT = " <<  runTime.deltaTValue() << endl;

// ************************************************************************* //
