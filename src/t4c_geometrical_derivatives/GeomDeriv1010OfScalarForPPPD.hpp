#ifndef GeomDeriv1010OfScalarForPPPD_hpp
#define GeomDeriv1010OfScalarForPPPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[PP|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_pppd: the integral geometrical derivatives buffer.
/// - Parameter buffer_spsd: the primitive integrals buffer.
/// - Parameter buffer_spdd: the primitive integrals buffer.
/// - Parameter buffer_dpsd: the primitive integrals buffer.
/// - Parameter buffer_dpdd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_pppd_0(CSimdArray<double>& buffer_1010_pppd,
                     const CSimdArray<double>& buffer_spsd,
                     const CSimdArray<double>& buffer_spdd,
                     const CSimdArray<double>& buffer_dpsd,
                     const CSimdArray<double>& buffer_dpdd,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForPPPD_hpp */
