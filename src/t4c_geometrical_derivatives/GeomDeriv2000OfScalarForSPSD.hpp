#ifndef GeomDeriv2000OfScalarForSPSD_hpp
#define GeomDeriv2000OfScalarForSPSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SP|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_spsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_spsd: the primitive integrals buffer.
/// - Parameter buffer_dpsd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_spsd_0(CSimdArray<double>& buffer_2000_spsd,
                     const CSimdArray<double>& buffer_spsd,
                     const CSimdArray<double>& buffer_dpsd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSPSD_hpp */
