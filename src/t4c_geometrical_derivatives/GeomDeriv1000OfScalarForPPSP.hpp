#ifndef GeomDeriv1000OfScalarForPPSP_hpp
#define GeomDeriv1000OfScalarForPPSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[PP|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_ppsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_spsp: the primitive integrals buffer.
/// - Parameter buffer_dpsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_ppsp_0(CSimdArray<double>& buffer_1000_ppsp,
                     const CSimdArray<double>& buffer_spsp,
                     const CSimdArray<double>& buffer_dpsp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForPPSP_hpp */
