#ifndef GeomDeriv1010OfScalarForPPPP_hpp
#define GeomDeriv1010OfScalarForPPPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[PP|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_pppp: the integral geometrical derivatives buffer.
/// - Parameter buffer_spsp: the primitive integrals buffer.
/// - Parameter buffer_spdp: the primitive integrals buffer.
/// - Parameter buffer_dpsp: the primitive integrals buffer.
/// - Parameter buffer_dpdp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_pppp_0(CSimdArray<double>& buffer_1010_pppp,
                     const CSimdArray<double>& buffer_spsp,
                     const CSimdArray<double>& buffer_spdp,
                     const CSimdArray<double>& buffer_dpsp,
                     const CSimdArray<double>& buffer_dpdp,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForPPPP_hpp */
