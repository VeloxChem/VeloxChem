#ifndef GeomDeriv1100OfScalarForSDSS_hpp
#define GeomDeriv1100OfScalarForSDSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[SD|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_sdss: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppss: the primitive integrals buffer.
/// - Parameter buffer_pfss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_sdss_0(CSimdArray<double>& buffer_1100_sdss,
                     const CSimdArray<double>& buffer_ppss,
                     const CSimdArray<double>& buffer_pfss,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForSDSS_hpp */
