#ifndef GeomDeriv2000OfScalarForSPPD_hpp
#define GeomDeriv2000OfScalarForSPPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SP|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_sppd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sppd: the primitive integrals buffer.
/// - Parameter buffer_dppd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_sppd_0(CSimdArray<double>& buffer_2000_sppd,
                     const CSimdArray<double>& buffer_sppd,
                     const CSimdArray<double>& buffer_dppd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSPPD_hpp */
