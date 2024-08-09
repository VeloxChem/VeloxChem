#ifndef GeomDeriv1100OfScalarForSSSP_hpp
#define GeomDeriv1100OfScalarForSSSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[SS|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_sssp: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_sssp_0(CSimdArray<double>& buffer_1100_sssp,
                     const CSimdArray<double>& buffer_ppsp,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForSSSP_hpp */
