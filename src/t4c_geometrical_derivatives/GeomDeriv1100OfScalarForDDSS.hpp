#ifndef GeomDeriv1100OfScalarForDDSS_hpp
#define GeomDeriv1100OfScalarForDDSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[DD|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_ddss: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppss: the primitive integrals buffer.
/// - Parameter buffer_pfss: the primitive integrals buffer.
/// - Parameter buffer_fpss: the primitive integrals buffer.
/// - Parameter buffer_ffss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_ddss_0(CSimdArray<double>& buffer_1100_ddss,
                     const CSimdArray<double>& buffer_ppss,
                     const CSimdArray<double>& buffer_pfss,
                     const CSimdArray<double>& buffer_fpss,
                     const CSimdArray<double>& buffer_ffss,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForDDSS_hpp */
