#ifndef GeomDeriv1100OfScalarForSSSS_hpp
#define GeomDeriv1100OfScalarForSSSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[SS|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_ssss: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_ssss_0(CSimdArray<double>& buffer_1100_ssss,
                     const CSimdArray<double>& buffer_ppss,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForSSSS_hpp */
