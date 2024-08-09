#ifndef GeomDeriv1100OfScalarForPDSD_hpp
#define GeomDeriv1100OfScalarForPDSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[PD|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_pdsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_spsd: the primitive integrals buffer.
/// - Parameter buffer_sfsd: the primitive integrals buffer.
/// - Parameter buffer_dpsd: the primitive integrals buffer.
/// - Parameter buffer_dfsd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_pdsd_0(CSimdArray<double>& buffer_1100_pdsd,
                     const CSimdArray<double>& buffer_spsd,
                     const CSimdArray<double>& buffer_sfsd,
                     const CSimdArray<double>& buffer_dpsd,
                     const CSimdArray<double>& buffer_dfsd,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForPDSD_hpp */
