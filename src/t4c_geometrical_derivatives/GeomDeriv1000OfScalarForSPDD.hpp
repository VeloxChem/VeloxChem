#ifndef GeomDeriv1000OfScalarForSPDD_hpp
#define GeomDeriv1000OfScalarForSPDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SP|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_spdd: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppdd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_spdd_0(CSimdArray<double>& buffer_1000_spdd,
                     const CSimdArray<double>& buffer_ppdd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSPDD_hpp */
