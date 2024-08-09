#ifndef GeomDeriv1100OfScalarForDDSD_hpp
#define GeomDeriv1100OfScalarForDDSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[DD|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_ddsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppsd: the primitive integrals buffer.
/// - Parameter buffer_pfsd: the primitive integrals buffer.
/// - Parameter buffer_fpsd: the primitive integrals buffer.
/// - Parameter buffer_ffsd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_ddsd_0(CSimdArray<double>& buffer_1100_ddsd,
                     const CSimdArray<double>& buffer_ppsd,
                     const CSimdArray<double>& buffer_pfsd,
                     const CSimdArray<double>& buffer_fpsd,
                     const CSimdArray<double>& buffer_ffsd,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForDDSD_hpp */
