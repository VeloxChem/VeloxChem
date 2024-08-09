#ifndef GeomDeriv1100OfScalarForSPSD_hpp
#define GeomDeriv1100OfScalarForSPSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[SP|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_spsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pssd: the primitive integrals buffer.
/// - Parameter buffer_pdsd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_spsd_0(CSimdArray<double>& buffer_1100_spsd,
                     const CSimdArray<double>& buffer_pssd,
                     const CSimdArray<double>& buffer_pdsd,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForSPSD_hpp */
