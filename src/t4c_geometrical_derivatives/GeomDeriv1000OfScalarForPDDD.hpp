#ifndef GeomDeriv1000OfScalarForPDDD_hpp
#define GeomDeriv1000OfScalarForPDDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[PD|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_pddd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sddd: the primitive integrals buffer.
/// - Parameter buffer_dddd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_pddd_0(CSimdArray<double>& buffer_1000_pddd,
                     const CSimdArray<double>& buffer_sddd,
                     const CSimdArray<double>& buffer_dddd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForPDDD_hpp */
