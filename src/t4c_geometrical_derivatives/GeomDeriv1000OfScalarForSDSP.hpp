#ifndef GeomDeriv1000OfScalarForSDSP_hpp
#define GeomDeriv1000OfScalarForSDSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SD|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_sdsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_sdsp_0(CSimdArray<double>& buffer_1000_sdsp,
                     const CSimdArray<double>& buffer_pdsp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSDSP_hpp */
