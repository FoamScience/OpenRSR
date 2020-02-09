#ifndef capillaryPressureModelBasesFwd_H
#define capillaryPressureModelBasesFwd_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<int Type> class capillaryPressureModelBase;

typedef capillaryPressureModelBase<2> base2PhasesPcModel; 
typedef capillaryPressureModelBase<3> base3PhasesPcModel;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
