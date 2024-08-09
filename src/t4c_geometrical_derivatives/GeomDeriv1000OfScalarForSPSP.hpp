#ifndef GeomDeriv1000OfScalarForSPSP_hpp
#define GeomDeriv1000OfScalarForSPSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SP|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_spsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_spsp_0(CSimdArray<double>& buffer_1000_spsp,
                     const CSimdArray<double>& buffer_ppsp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSPSP_hpp */
