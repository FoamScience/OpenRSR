{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(-Mf, p) + fvc::div(phiG)
        // Capillary flux
        + fvc::div(phiPc)
        // Wells contribution
        + wModel->source(canPhase.name())
        + wModel->source(noncanPhase.name())
        //+ SrcProd -SrcInj
        //+ ewSource + eoSource
        //+ fvm::Sp(iwSource, p) + fvm::Sp(ioSource, p)
    );

    pEqn.solve();

    phiP = pEqn.flux();
    //Equivqlent to: phiP = -Mf * fvc::snGrad(p) * mesh.magSf();

}
