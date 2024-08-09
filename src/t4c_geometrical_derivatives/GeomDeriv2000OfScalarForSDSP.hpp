#ifndef GeomDeriv2000OfScalarForSDSP_hpp
#define GeomDeriv2000OfScalarForSDSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SD|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_sdsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdsp: the primitive integrals buffer.
/// - Parameter buffer_ddsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_sdsp_0(CSimdArray<double>& buffer_2000_sdsp,
                     const CSimdArray<double>& buffer_sdsp,
                     const CSimdArray<double>& buffer_ddsp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSDSP_hpp */
