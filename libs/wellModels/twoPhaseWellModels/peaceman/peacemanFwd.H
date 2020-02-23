#ifndef peacemanFwd_H
#define peacemanFwd_H

#include "parameterTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class KType, class MuRhoType> class peaceman;

typedef peaceman<Iso, Incompressible>    derivedIsoIncompressiblePeaceman;
typedef peaceman<Iso, Compressible>      derivedIsoCompressiblePeaceman;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
