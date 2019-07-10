{
    fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mnf, p) + fvc::div(phinG)
            // Capillary flux
            //+ fvc::div(phiPc)
            // Wells contribution
            //+ SrcProd -SrcInj
            + eoSource //+ ewSource
            + fvm::Sp(ioSource, p) //+ fvm::Sp(iwSource, p)
        );

    fvScalarMatrix SwInp
        (
            -fvm::ddt(oStorage, phasew.alpha())
        );

    resEqn.insertEquation(1,pEqn);
    resEqn.insertEquationCoupling(1,0,SwInp);
    

    //phiP = pEqn.flux();
    //eq. to: phiP = -Mf*fvc::snGrad(p) * mesh.magSf();

}
