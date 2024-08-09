#ifndef GeomDeriv1100OfScalarForPPPD_hpp
#define GeomDeriv1100OfScalarForPPPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[PP|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_pppd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sspd: the primitive integrals buffer.
/// - Parameter buffer_sdpd: the primitive integrals buffer.
/// - Parameter buffer_dspd: the primitive integrals buffer.
/// - Parameter buffer_ddpd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_pppd_0(CSimdArray<double>& buffer_1100_pppd,
                     const CSimdArray<double>& buffer_sspd,
                     const CSimdArray<double>& buffer_sdpd,
                     const CSimdArray<double>& buffer_dspd,
                     const CSimdArray<double>& buffer_ddpd,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForPPPD_hpp */
