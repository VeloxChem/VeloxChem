#ifndef GeomDeriv2000OfScalarForPDDD_hpp
#define GeomDeriv2000OfScalarForPDDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[PD|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_pddd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pddd: the primitive integrals buffer.
/// - Parameter buffer_fddd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_pddd_0(CSimdArray<double>& buffer_2000_pddd,
                     const CSimdArray<double>& buffer_pddd,
                     const CSimdArray<double>& buffer_fddd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForPDDD_hpp */
