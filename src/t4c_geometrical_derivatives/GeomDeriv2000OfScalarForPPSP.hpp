#ifndef GeomDeriv2000OfScalarForPPSP_hpp
#define GeomDeriv2000OfScalarForPPSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[PP|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_ppsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppsp: the primitive integrals buffer.
/// - Parameter buffer_fpsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_ppsp_0(CSimdArray<double>& buffer_2000_ppsp,
                     const CSimdArray<double>& buffer_ppsp,
                     const CSimdArray<double>& buffer_fpsp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForPPSP_hpp */
