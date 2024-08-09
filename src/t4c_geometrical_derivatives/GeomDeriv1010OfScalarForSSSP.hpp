#ifndef GeomDeriv1010OfScalarForSSSP_hpp
#define GeomDeriv1010OfScalarForSSSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[SS|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_sssp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pspp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_sssp_0(CSimdArray<double>& buffer_1010_sssp,
                     const CSimdArray<double>& buffer_pspp,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForSSSP_hpp */
