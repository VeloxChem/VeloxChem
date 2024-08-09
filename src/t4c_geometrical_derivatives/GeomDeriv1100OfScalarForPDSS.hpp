#ifndef GeomDeriv1100OfScalarForPDSS_hpp
#define GeomDeriv1100OfScalarForPDSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[PD|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_pdss: the integral geometrical derivatives buffer.
/// - Parameter buffer_spss: the primitive integrals buffer.
/// - Parameter buffer_sfss: the primitive integrals buffer.
/// - Parameter buffer_dpss: the primitive integrals buffer.
/// - Parameter buffer_dfss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_pdss_0(CSimdArray<double>& buffer_1100_pdss,
                     const CSimdArray<double>& buffer_spss,
                     const CSimdArray<double>& buffer_sfss,
                     const CSimdArray<double>& buffer_dpss,
                     const CSimdArray<double>& buffer_dfss,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForPDSS_hpp */
