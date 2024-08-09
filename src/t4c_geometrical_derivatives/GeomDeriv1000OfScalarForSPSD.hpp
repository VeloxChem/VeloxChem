#ifndef GeomDeriv1000OfScalarForSPSD_hpp
#define GeomDeriv1000OfScalarForSPSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SP|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_spsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppsd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_spsd_0(CSimdArray<double>& buffer_1000_spsd,
                     const CSimdArray<double>& buffer_ppsd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSPSD_hpp */
