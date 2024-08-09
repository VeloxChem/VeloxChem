#ifndef GeomDeriv1100OfScalarForPPDD_hpp
#define GeomDeriv1100OfScalarForPPDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[PP|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_ppdd: the integral geometrical derivatives buffer.
/// - Parameter buffer_ssdd: the primitive integrals buffer.
/// - Parameter buffer_sddd: the primitive integrals buffer.
/// - Parameter buffer_dsdd: the primitive integrals buffer.
/// - Parameter buffer_dddd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_ppdd_0(CSimdArray<double>& buffer_1100_ppdd,
                     const CSimdArray<double>& buffer_ssdd,
                     const CSimdArray<double>& buffer_sddd,
                     const CSimdArray<double>& buffer_dsdd,
                     const CSimdArray<double>& buffer_dddd,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForPPDD_hpp */
