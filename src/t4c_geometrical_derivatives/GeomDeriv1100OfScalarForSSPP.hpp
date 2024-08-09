#ifndef GeomDeriv1100OfScalarForSSPP_hpp
#define GeomDeriv1100OfScalarForSSPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[SS|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_sspp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pppp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_sspp_0(CSimdArray<double>& buffer_1100_sspp,
                     const CSimdArray<double>& buffer_pppp,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForSSPP_hpp */
