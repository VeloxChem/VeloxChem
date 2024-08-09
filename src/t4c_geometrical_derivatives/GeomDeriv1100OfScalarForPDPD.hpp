#ifndef GeomDeriv1100OfScalarForPDPD_hpp
#define GeomDeriv1100OfScalarForPDPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[PD|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_pdpd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sppd: the primitive integrals buffer.
/// - Parameter buffer_sfpd: the primitive integrals buffer.
/// - Parameter buffer_dppd: the primitive integrals buffer.
/// - Parameter buffer_dfpd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_pdpd_0(CSimdArray<double>& buffer_1100_pdpd,
                     const CSimdArray<double>& buffer_sppd,
                     const CSimdArray<double>& buffer_sfpd,
                     const CSimdArray<double>& buffer_dppd,
                     const CSimdArray<double>& buffer_dfpd,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForPDPD_hpp */
