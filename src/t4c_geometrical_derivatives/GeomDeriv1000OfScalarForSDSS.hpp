#ifndef GeomDeriv1000OfScalarForSDSS_hpp
#define GeomDeriv1000OfScalarForSDSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SD|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_sdss: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_sdss_0(CSimdArray<double>& buffer_1000_sdss,
                     const CSimdArray<double>& buffer_pdss,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSDSS_hpp */
