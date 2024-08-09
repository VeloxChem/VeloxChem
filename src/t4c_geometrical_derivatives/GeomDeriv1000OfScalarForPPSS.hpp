#ifndef GeomDeriv1000OfScalarForPPSS_hpp
#define GeomDeriv1000OfScalarForPPSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[PP|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_ppss: the integral geometrical derivatives buffer.
/// - Parameter buffer_spss: the primitive integrals buffer.
/// - Parameter buffer_dpss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_ppss_0(CSimdArray<double>& buffer_1000_ppss,
                     const CSimdArray<double>& buffer_spss,
                     const CSimdArray<double>& buffer_dpss,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForPPSS_hpp */
