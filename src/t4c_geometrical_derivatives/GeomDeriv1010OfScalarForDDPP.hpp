#ifndef GeomDeriv1010OfScalarForDDPP_hpp
#define GeomDeriv1010OfScalarForDDPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[DD|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_ddpp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdsp: the primitive integrals buffer.
/// - Parameter buffer_pddp: the primitive integrals buffer.
/// - Parameter buffer_fdsp: the primitive integrals buffer.
/// - Parameter buffer_fddp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_ddpp_0(CSimdArray<double>& buffer_1010_ddpp,
                     const CSimdArray<double>& buffer_pdsp,
                     const CSimdArray<double>& buffer_pddp,
                     const CSimdArray<double>& buffer_fdsp,
                     const CSimdArray<double>& buffer_fddp,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForDDPP_hpp */
