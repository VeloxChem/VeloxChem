#ifndef GeomDeriv1100OfScalarForPPSS_hpp
#define GeomDeriv1100OfScalarForPPSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[PP|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_ppss: the integral geometrical derivatives buffer.
/// - Parameter buffer_ssss: the primitive integrals buffer.
/// - Parameter buffer_sdss: the primitive integrals buffer.
/// - Parameter buffer_dsss: the primitive integrals buffer.
/// - Parameter buffer_ddss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_ppss_0(CSimdArray<double>& buffer_1100_ppss,
                     const CSimdArray<double>& buffer_ssss,
                     const CSimdArray<double>& buffer_sdss,
                     const CSimdArray<double>& buffer_dsss,
                     const CSimdArray<double>& buffer_ddss,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForPPSS_hpp */
