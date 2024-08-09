#ifndef GeomDeriv1100OfScalarForPPPP_hpp
#define GeomDeriv1100OfScalarForPPPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[PP|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_pppp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sspp: the primitive integrals buffer.
/// - Parameter buffer_sdpp: the primitive integrals buffer.
/// - Parameter buffer_dspp: the primitive integrals buffer.
/// - Parameter buffer_ddpp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_pppp_0(CSimdArray<double>& buffer_1100_pppp,
                     const CSimdArray<double>& buffer_sspp,
                     const CSimdArray<double>& buffer_sdpp,
                     const CSimdArray<double>& buffer_dspp,
                     const CSimdArray<double>& buffer_ddpp,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForPPPP_hpp */
