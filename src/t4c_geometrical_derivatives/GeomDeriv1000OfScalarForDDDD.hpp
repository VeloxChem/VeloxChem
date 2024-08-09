#ifndef GeomDeriv1000OfScalarForDDDD_hpp
#define GeomDeriv1000OfScalarForDDDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[DD|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_dddd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pddd: the primitive integrals buffer.
/// - Parameter buffer_fddd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_dddd_0(CSimdArray<double>& buffer_1000_dddd,
                     const CSimdArray<double>& buffer_pddd,
                     const CSimdArray<double>& buffer_fddd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForDDDD_hpp */
