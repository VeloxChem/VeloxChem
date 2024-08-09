#ifndef GeomDeriv1000OfScalarForSDDD_hpp
#define GeomDeriv1000OfScalarForSDDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SD|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_sddd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pddd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_sddd_0(CSimdArray<double>& buffer_1000_sddd,
                     const CSimdArray<double>& buffer_pddd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSDDD_hpp */
