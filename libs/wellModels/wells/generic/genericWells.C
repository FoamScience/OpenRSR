
#include "generic.C"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace wells
{
    // Isotropic and anisotropic spesializations 
    typedef generic<volScalarField> isoGeneric;
    //typedef generic<volTensorField> anisoGeneric;

    // Add both versions to RTS table
    defineTemplateTypeNameAndDebug(isoGeneric, 0);
    //defineTemplateTypeNameAndDebug(anisoGeneric, 0);
    
    addToRunTimeSelectionTable
    (
        well,
        isoGeneric,
        dictionary
    );
    //addToRunTimeSelectionTable
    //(
    //    well,
    //    anisoGeneric,
    //    dictionary
    //);
}
}
