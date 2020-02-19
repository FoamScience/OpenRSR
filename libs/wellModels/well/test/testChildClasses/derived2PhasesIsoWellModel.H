// A child well model for isotropic 2-phase wells
// Providing a sample operator()
struct wellModelDerived : public base2IsoWellModel
{
    TypeName("wellModelDerived");
    wellModelDerived(
            const word& name,
            const dictionary& wellProperties,
            const fvMesh& mesh
            ): base2IsoWellModel(name, wellProperties, mesh) {}
    virtual void correct(){}
    virtual void operator()(const word& wellName) const
    {
        if (foundObject<base2IsoWell>(wellName))
        {
            const base2IsoWell& well = lookupObject<base2IsoWell>(wellName);
            // Operate on well data members
            well.setEquivalentRadius(15);
        }
    }
};

#include "addToRunTimeSelectionTable.H"
defineTypeNameAndDebug(wellModelDerived, 0);
addToRunTimeSelectionTable(base2IsoWellModel, wellModelDerived, dictionary);
