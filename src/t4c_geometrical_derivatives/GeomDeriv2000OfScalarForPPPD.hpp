#ifndef GeomDeriv2000OfScalarForPPPD_hpp
#define GeomDeriv2000OfScalarForPPPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[PP|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_pppd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pppd: the primitive integrals buffer.
/// - Parameter buffer_fppd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_pppd_0(CSimdArray<double>& buffer_2000_pppd,
                     const CSimdArray<double>& buffer_pppd,
                     const CSimdArray<double>& buffer_fppd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForPPPD_hpp */
