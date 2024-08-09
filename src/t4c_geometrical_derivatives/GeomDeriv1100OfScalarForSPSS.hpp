#ifndef GeomDeriv1100OfScalarForSPSS_hpp
#define GeomDeriv1100OfScalarForSPSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[SP|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_spss: the integral geometrical derivatives buffer.
/// - Parameter buffer_psss: the primitive integrals buffer.
/// - Parameter buffer_pdss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_spss_0(CSimdArray<double>& buffer_1100_spss,
                     const CSimdArray<double>& buffer_psss,
                     const CSimdArray<double>& buffer_pdss,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForSPSS_hpp */
