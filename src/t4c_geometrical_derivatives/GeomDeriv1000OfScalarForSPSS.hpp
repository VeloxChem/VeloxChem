#ifndef GeomDeriv1000OfScalarForSPSS_hpp
#define GeomDeriv1000OfScalarForSPSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SP|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_spss: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_spss_0(CSimdArray<double>& buffer_1000_spss,
                     const CSimdArray<double>& buffer_ppss,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSPSS_hpp */
