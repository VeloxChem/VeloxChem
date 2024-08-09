#ifndef GeomDeriv2000OfScalarForDDDD_hpp
#define GeomDeriv2000OfScalarForDDDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[DD|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_dddd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sddd: the primitive integrals buffer.
/// - Parameter buffer_dddd: the primitive integrals buffer.
/// - Parameter buffer_gddd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_dddd_0(CSimdArray<double>& buffer_2000_dddd,
                     const CSimdArray<double>& buffer_sddd,
                     const CSimdArray<double>& buffer_dddd,
                     const CSimdArray<double>& buffer_gddd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForDDDD_hpp */
