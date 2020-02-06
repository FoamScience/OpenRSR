#ifndef relativePermeabilityModelBasesFwd_H
#define relativePermeabilityModelBasesFwd_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<int Type> class relativePermeabilityModelBase;

typedef relativePermeabilityModelBase<2> base2PhasesKrModel; 
typedef relativePermeabilityModelBase<3> base3PhasesKrModel;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
