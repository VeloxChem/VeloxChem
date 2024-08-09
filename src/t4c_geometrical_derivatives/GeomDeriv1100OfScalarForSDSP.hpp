#ifndef GeomDeriv1100OfScalarForSDSP_hpp
#define GeomDeriv1100OfScalarForSDSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[SD|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_sdsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppsp: the primitive integrals buffer.
/// - Parameter buffer_pfsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_sdsp_0(CSimdArray<double>& buffer_1100_sdsp,
                     const CSimdArray<double>& buffer_ppsp,
                     const CSimdArray<double>& buffer_pfsp,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForSDSP_hpp */
