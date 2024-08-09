#ifndef GeomDeriv1100OfScalarForPDPP_hpp
#define GeomDeriv1100OfScalarForPDPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[PD|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_pdpp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sppp: the primitive integrals buffer.
/// - Parameter buffer_sfpp: the primitive integrals buffer.
/// - Parameter buffer_dppp: the primitive integrals buffer.
/// - Parameter buffer_dfpp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_pdpp_0(CSimdArray<double>& buffer_1100_pdpp,
                     const CSimdArray<double>& buffer_sppp,
                     const CSimdArray<double>& buffer_sfpp,
                     const CSimdArray<double>& buffer_dppp,
                     const CSimdArray<double>& buffer_dfpp,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForPDPP_hpp */
