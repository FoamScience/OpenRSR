{
    fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mf, p) + fvc::div(phiG)
            // Capillary flux
            + fvc::div(phiPc)
            // Wells contribution
            + wModel->source(canPhase.name())
            + wModel->source(noncanPhase.name())
        );

    pEqn.solve();

    phiP = pEqn.flux();
    //phiP = -Mf*fvc::snGrad(p) * mesh.magSf();
}
