#ifndef GeomDeriv1100OfScalarForPDDD_hpp
#define GeomDeriv1100OfScalarForPDDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[PD|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_pddd: the integral geometrical derivatives buffer.
/// - Parameter buffer_spdd: the primitive integrals buffer.
/// - Parameter buffer_sfdd: the primitive integrals buffer.
/// - Parameter buffer_dpdd: the primitive integrals buffer.
/// - Parameter buffer_dfdd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_pddd_0(CSimdArray<double>& buffer_1100_pddd,
                     const CSimdArray<double>& buffer_spdd,
                     const CSimdArray<double>& buffer_sfdd,
                     const CSimdArray<double>& buffer_dpdd,
                     const CSimdArray<double>& buffer_dfdd,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForPDDD_hpp */
