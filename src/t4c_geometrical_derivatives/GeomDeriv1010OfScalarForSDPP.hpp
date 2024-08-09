#ifndef GeomDeriv1010OfScalarForSDPP_hpp
#define GeomDeriv1010OfScalarForSDPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[SD|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_sdpp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdsp: the primitive integrals buffer.
/// - Parameter buffer_pddp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_sdpp_0(CSimdArray<double>& buffer_1010_sdpp,
                     const CSimdArray<double>& buffer_pdsp,
                     const CSimdArray<double>& buffer_pddp,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForSDPP_hpp */
