#ifndef GeomDeriv2000OfScalarForDDPP_hpp
#define GeomDeriv2000OfScalarForDDPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[DD|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_ddpp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdpp: the primitive integrals buffer.
/// - Parameter buffer_ddpp: the primitive integrals buffer.
/// - Parameter buffer_gdpp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_ddpp_0(CSimdArray<double>& buffer_2000_ddpp,
                     const CSimdArray<double>& buffer_sdpp,
                     const CSimdArray<double>& buffer_ddpp,
                     const CSimdArray<double>& buffer_gdpp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForDDPP_hpp */
