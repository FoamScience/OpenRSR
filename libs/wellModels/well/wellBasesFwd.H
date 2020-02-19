#ifndef wellBasesFwd_H
#define wellBasesFwd_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class KType,int nPhases> class wellBase;
template<class KType,int nPhases> class wellModelBase;

typedef wellBase<Iso,2>   base2IsoWell;
typedef wellBase<Aniso,2> base2AnisoWell;

typedef wellModelBase<Iso,2>   base2IsoWellModel;
typedef wellModelBase<Aniso,2> base2AnisoWellModel;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
