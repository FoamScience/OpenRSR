Rw = mag
    ( 
        fvc::ddt(porosity, phasew.alpha()) 
        + fvc::div(phiwG) + fvc::div(phiPc)
        + fvc::laplacian(-Mwf,p)
        + iwSource*p
        + ewSource
    );
Rn = mag
    (
            fvc::laplacian(-Mnf, p) + fvc::div(phinG)
            // Capillary flux
            //+ fvc::div(phiPc)
            // Wells contribution
            //+ SrcProd -SrcInj
            + eoSource //+ ewSource
            + ioSource*p//+ fvm::Sp(iwSource, p)
            - fvc::ddt(oStorage, phasew.alpha())
    );
nRw = Rw/porosity.internalField()/mesh.V().field()/rFVFw.internalField();
nRn = Rn/porosity.internalField()/mesh.V().field()/rFVFn.internalField();
CNVw = runTime.deltaT().value()*gMax(nRw);
CNVn = runTime.deltaT().value()*gMax(nRn);

Info<< "CNV numbers: " << phasew.name() << ".CNV " << CNVw 
    << ", " << phasen.name() << ".CNV " << CNVn << endl;
