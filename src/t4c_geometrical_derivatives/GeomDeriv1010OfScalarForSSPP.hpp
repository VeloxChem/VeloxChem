#ifndef GeomDeriv1010OfScalarForSSPP_hpp
#define GeomDeriv1010OfScalarForSSPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[SS|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_sspp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pssp: the primitive integrals buffer.
/// - Parameter buffer_psdp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_sspp_0(CSimdArray<double>& buffer_1010_sspp,
                     const CSimdArray<double>& buffer_pssp,
                     const CSimdArray<double>& buffer_psdp,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForSSPP_hpp */
