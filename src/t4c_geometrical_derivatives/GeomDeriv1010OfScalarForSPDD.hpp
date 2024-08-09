#ifndef GeomDeriv1010OfScalarForSPDD_hpp
#define GeomDeriv1010OfScalarForSPDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[SP|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_spdd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pppd: the primitive integrals buffer.
/// - Parameter buffer_ppfd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_spdd_0(CSimdArray<double>& buffer_1010_spdd,
                     const CSimdArray<double>& buffer_pppd,
                     const CSimdArray<double>& buffer_ppfd,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForSPDD_hpp */
