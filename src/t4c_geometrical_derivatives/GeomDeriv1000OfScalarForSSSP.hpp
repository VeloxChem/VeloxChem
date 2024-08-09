#ifndef GeomDeriv1000OfScalarForSSSP_hpp
#define GeomDeriv1000OfScalarForSSSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SS|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_sssp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pssp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_sssp_0(CSimdArray<double>& buffer_1000_sssp,
                     const CSimdArray<double>& buffer_pssp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSSSP_hpp */
