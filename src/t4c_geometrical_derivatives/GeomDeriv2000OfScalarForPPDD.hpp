#ifndef GeomDeriv2000OfScalarForPPDD_hpp
#define GeomDeriv2000OfScalarForPPDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[PP|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_ppdd: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppdd: the primitive integrals buffer.
/// - Parameter buffer_fpdd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_ppdd_0(CSimdArray<double>& buffer_2000_ppdd,
                     const CSimdArray<double>& buffer_ppdd,
                     const CSimdArray<double>& buffer_fpdd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForPPDD_hpp */
