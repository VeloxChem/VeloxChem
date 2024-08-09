#ifndef GeomDeriv2000OfScalarForSDPP_hpp
#define GeomDeriv2000OfScalarForSDPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SD|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_sdpp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdpp: the primitive integrals buffer.
/// - Parameter buffer_ddpp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_sdpp_0(CSimdArray<double>& buffer_2000_sdpp,
                     const CSimdArray<double>& buffer_sdpp,
                     const CSimdArray<double>& buffer_ddpp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSDPP_hpp */
