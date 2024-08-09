#ifndef GeomDeriv1010OfScalarForSPSP_hpp
#define GeomDeriv1010OfScalarForSPSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[SP|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_spsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pppp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_spsp_0(CSimdArray<double>& buffer_1010_spsp,
                     const CSimdArray<double>& buffer_pppp,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForSPSP_hpp */
