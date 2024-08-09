#ifndef GeomDeriv1010OfScalarForSDSP_hpp
#define GeomDeriv1010OfScalarForSDSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[SD|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_sdsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdpp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_sdsp_0(CSimdArray<double>& buffer_1010_sdsp,
                     const CSimdArray<double>& buffer_pdpp,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForSDSP_hpp */
