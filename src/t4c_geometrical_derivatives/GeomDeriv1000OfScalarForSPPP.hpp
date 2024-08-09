#ifndef GeomDeriv1000OfScalarForSPPP_hpp
#define GeomDeriv1000OfScalarForSPPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SP|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_sppp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pppp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_sppp_0(CSimdArray<double>& buffer_1000_sppp,
                     const CSimdArray<double>& buffer_pppp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSPPP_hpp */
