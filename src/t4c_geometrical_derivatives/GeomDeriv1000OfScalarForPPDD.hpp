#ifndef GeomDeriv1000OfScalarForPPDD_hpp
#define GeomDeriv1000OfScalarForPPDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[PP|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_ppdd: the integral geometrical derivatives buffer.
/// - Parameter buffer_spdd: the primitive integrals buffer.
/// - Parameter buffer_dpdd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_ppdd_0(CSimdArray<double>& buffer_1000_ppdd,
                     const CSimdArray<double>& buffer_spdd,
                     const CSimdArray<double>& buffer_dpdd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForPPDD_hpp */
