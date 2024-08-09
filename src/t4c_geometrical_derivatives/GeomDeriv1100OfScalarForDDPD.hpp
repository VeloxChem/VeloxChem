#ifndef GeomDeriv1100OfScalarForDDPD_hpp
#define GeomDeriv1100OfScalarForDDPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[DD|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_ddpd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pppd: the primitive integrals buffer.
/// - Parameter buffer_pfpd: the primitive integrals buffer.
/// - Parameter buffer_fppd: the primitive integrals buffer.
/// - Parameter buffer_ffpd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_ddpd_0(CSimdArray<double>& buffer_1100_ddpd,
                     const CSimdArray<double>& buffer_pppd,
                     const CSimdArray<double>& buffer_pfpd,
                     const CSimdArray<double>& buffer_fppd,
                     const CSimdArray<double>& buffer_ffpd,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForDDPD_hpp */
