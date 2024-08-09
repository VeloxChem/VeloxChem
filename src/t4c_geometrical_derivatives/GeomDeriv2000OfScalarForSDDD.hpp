#ifndef GeomDeriv2000OfScalarForSDDD_hpp
#define GeomDeriv2000OfScalarForSDDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SD|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_sddd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sddd: the primitive integrals buffer.
/// - Parameter buffer_dddd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_sddd_0(CSimdArray<double>& buffer_2000_sddd,
                     const CSimdArray<double>& buffer_sddd,
                     const CSimdArray<double>& buffer_dddd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSDDD_hpp */
