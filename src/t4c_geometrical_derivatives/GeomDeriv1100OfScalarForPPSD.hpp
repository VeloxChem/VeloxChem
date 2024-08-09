#ifndef GeomDeriv1100OfScalarForPPSD_hpp
#define GeomDeriv1100OfScalarForPPSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[PP|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_ppsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sssd: the primitive integrals buffer.
/// - Parameter buffer_sdsd: the primitive integrals buffer.
/// - Parameter buffer_dssd: the primitive integrals buffer.
/// - Parameter buffer_ddsd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_ppsd_0(CSimdArray<double>& buffer_1100_ppsd,
                     const CSimdArray<double>& buffer_sssd,
                     const CSimdArray<double>& buffer_sdsd,
                     const CSimdArray<double>& buffer_dssd,
                     const CSimdArray<double>& buffer_ddsd,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForPPSD_hpp */
