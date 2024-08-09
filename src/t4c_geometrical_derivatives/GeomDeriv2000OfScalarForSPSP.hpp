#ifndef GeomDeriv2000OfScalarForSPSP_hpp
#define GeomDeriv2000OfScalarForSPSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SP|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_spsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_spsp: the primitive integrals buffer.
/// - Parameter buffer_dpsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_spsp_0(CSimdArray<double>& buffer_2000_spsp,
                     const CSimdArray<double>& buffer_spsp,
                     const CSimdArray<double>& buffer_dpsp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSPSP_hpp */
