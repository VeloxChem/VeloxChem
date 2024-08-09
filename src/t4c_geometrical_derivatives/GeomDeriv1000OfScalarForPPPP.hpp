#ifndef GeomDeriv1000OfScalarForPPPP_hpp
#define GeomDeriv1000OfScalarForPPPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[PP|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_pppp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sppp: the primitive integrals buffer.
/// - Parameter buffer_dppp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_pppp_0(CSimdArray<double>& buffer_1000_pppp,
                     const CSimdArray<double>& buffer_sppp,
                     const CSimdArray<double>& buffer_dppp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForPPPP_hpp */
