#ifndef GeomDeriv1100OfScalarForPPSP_hpp
#define GeomDeriv1100OfScalarForPPSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[PP|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_ppsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sssp: the primitive integrals buffer.
/// - Parameter buffer_sdsp: the primitive integrals buffer.
/// - Parameter buffer_dssp: the primitive integrals buffer.
/// - Parameter buffer_ddsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_ppsp_0(CSimdArray<double>& buffer_1100_ppsp,
                     const CSimdArray<double>& buffer_sssp,
                     const CSimdArray<double>& buffer_sdsp,
                     const CSimdArray<double>& buffer_dssp,
                     const CSimdArray<double>& buffer_ddsp,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForPPSP_hpp */
