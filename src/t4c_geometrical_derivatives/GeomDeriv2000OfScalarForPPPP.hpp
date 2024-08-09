#ifndef GeomDeriv2000OfScalarForPPPP_hpp
#define GeomDeriv2000OfScalarForPPPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[PP|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_pppp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pppp: the primitive integrals buffer.
/// - Parameter buffer_fppp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_pppp_0(CSimdArray<double>& buffer_2000_pppp,
                     const CSimdArray<double>& buffer_pppp,
                     const CSimdArray<double>& buffer_fppp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForPPPP_hpp */
