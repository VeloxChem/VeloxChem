#ifndef GeomDeriv1100OfScalarForSPSP_hpp
#define GeomDeriv1100OfScalarForSPSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[SP|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_spsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pssp: the primitive integrals buffer.
/// - Parameter buffer_pdsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_spsp_0(CSimdArray<double>& buffer_1100_spsp,
                     const CSimdArray<double>& buffer_pssp,
                     const CSimdArray<double>& buffer_pdsp,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForSPSP_hpp */
