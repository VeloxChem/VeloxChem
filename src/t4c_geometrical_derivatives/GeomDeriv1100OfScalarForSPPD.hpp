#ifndef GeomDeriv1100OfScalarForSPPD_hpp
#define GeomDeriv1100OfScalarForSPPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[SP|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_sppd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pspd: the primitive integrals buffer.
/// - Parameter buffer_pdpd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_sppd_0(CSimdArray<double>& buffer_1100_sppd,
                     const CSimdArray<double>& buffer_pspd,
                     const CSimdArray<double>& buffer_pdpd,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForSPPD_hpp */
