#ifndef GeomDeriv1100OfScalarForSPPP_hpp
#define GeomDeriv1100OfScalarForSPPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[SP|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_sppp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pspp: the primitive integrals buffer.
/// - Parameter buffer_pdpp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_sppp_0(CSimdArray<double>& buffer_1100_sppp,
                     const CSimdArray<double>& buffer_pspp,
                     const CSimdArray<double>& buffer_pdpp,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForSPPP_hpp */
