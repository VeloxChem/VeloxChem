#ifndef GeomDeriv2000OfScalarForSDSS_hpp
#define GeomDeriv2000OfScalarForSDSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SD|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_sdss: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdss: the primitive integrals buffer.
/// - Parameter buffer_ddss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_sdss_0(CSimdArray<double>& buffer_2000_sdss,
                     const CSimdArray<double>& buffer_sdss,
                     const CSimdArray<double>& buffer_ddss,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSDSS_hpp */
