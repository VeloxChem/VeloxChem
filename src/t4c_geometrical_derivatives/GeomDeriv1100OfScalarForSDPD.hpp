#ifndef GeomDeriv1100OfScalarForSDPD_hpp
#define GeomDeriv1100OfScalarForSDPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[SD|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_sdpd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pppd: the primitive integrals buffer.
/// - Parameter buffer_pfpd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_sdpd_0(CSimdArray<double>& buffer_1100_sdpd,
                     const CSimdArray<double>& buffer_pppd,
                     const CSimdArray<double>& buffer_pfpd,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForSDPD_hpp */
