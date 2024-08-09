#ifndef GeomDeriv1100OfScalarForSSPD_hpp
#define GeomDeriv1100OfScalarForSSPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[SS|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_sspd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pppd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_sspd_0(CSimdArray<double>& buffer_1100_sspd,
                     const CSimdArray<double>& buffer_pppd,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForSSPD_hpp */
