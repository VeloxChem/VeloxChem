#ifndef GeomDeriv1100OfScalarForSDPP_hpp
#define GeomDeriv1100OfScalarForSDPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[SD|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_sdpp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pppp: the primitive integrals buffer.
/// - Parameter buffer_pfpp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_sdpp_0(CSimdArray<double>& buffer_1100_sdpp,
                     const CSimdArray<double>& buffer_pppp,
                     const CSimdArray<double>& buffer_pfpp,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForSDPP_hpp */
