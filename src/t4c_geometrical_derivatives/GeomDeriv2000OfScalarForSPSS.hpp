#ifndef GeomDeriv2000OfScalarForSPSS_hpp
#define GeomDeriv2000OfScalarForSPSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SP|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_spss: the integral geometrical derivatives buffer.
/// - Parameter buffer_spss: the primitive integrals buffer.
/// - Parameter buffer_dpss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_spss_0(CSimdArray<double>& buffer_2000_spss,
                     const CSimdArray<double>& buffer_spss,
                     const CSimdArray<double>& buffer_dpss,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSPSS_hpp */
