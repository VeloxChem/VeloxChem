#ifndef GeomDeriv2000OfScalarForSPDD_hpp
#define GeomDeriv2000OfScalarForSPDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SP|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_spdd: the integral geometrical derivatives buffer.
/// - Parameter buffer_spdd: the primitive integrals buffer.
/// - Parameter buffer_dpdd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_spdd_0(CSimdArray<double>& buffer_2000_spdd,
                     const CSimdArray<double>& buffer_spdd,
                     const CSimdArray<double>& buffer_dpdd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSPDD_hpp */
