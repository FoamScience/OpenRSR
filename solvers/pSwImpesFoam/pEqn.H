{
    fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mf, p) + fvc::div(phiG)
            // Capillary flux
            + fvc::div(phiPc)
            // Wells contribution
            //+ SrcProd -SrcInj
            + ewSource + eoSource
            + fvm::Sp(iwSource, p) + fvm::Sp(ioSource, p)
        );

    pEqn.solve();

    phiP = pEqn.flux();
    //eq. to: phiP = -Mf*fvc::snGrad(p) * mesh.magSf();

}
