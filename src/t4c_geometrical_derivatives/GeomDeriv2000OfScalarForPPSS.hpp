#ifndef GeomDeriv2000OfScalarForPPSS_hpp
#define GeomDeriv2000OfScalarForPPSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[PP|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_ppss: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppss: the primitive integrals buffer.
/// - Parameter buffer_fpss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_ppss_0(CSimdArray<double>& buffer_2000_ppss,
                     const CSimdArray<double>& buffer_ppss,
                     const CSimdArray<double>& buffer_fpss,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForPPSS_hpp */
