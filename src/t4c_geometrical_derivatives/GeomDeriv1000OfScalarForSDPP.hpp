#ifndef GeomDeriv1000OfScalarForSDPP_hpp
#define GeomDeriv1000OfScalarForSDPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SD|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_sdpp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdpp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_sdpp_0(CSimdArray<double>& buffer_1000_sdpp,
                     const CSimdArray<double>& buffer_pdpp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSDPP_hpp */
