#ifndef GeomDeriv2000OfScalarForSPPP_hpp
#define GeomDeriv2000OfScalarForSPPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SP|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_sppp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sppp: the primitive integrals buffer.
/// - Parameter buffer_dppp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_sppp_0(CSimdArray<double>& buffer_2000_sppp,
                     const CSimdArray<double>& buffer_sppp,
                     const CSimdArray<double>& buffer_dppp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSPPP_hpp */
