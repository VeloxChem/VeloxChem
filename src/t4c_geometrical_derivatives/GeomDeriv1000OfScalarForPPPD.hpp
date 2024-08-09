#ifndef GeomDeriv1000OfScalarForPPPD_hpp
#define GeomDeriv1000OfScalarForPPPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[PP|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_pppd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sppd: the primitive integrals buffer.
/// - Parameter buffer_dppd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_pppd_0(CSimdArray<double>& buffer_1000_pppd,
                     const CSimdArray<double>& buffer_sppd,
                     const CSimdArray<double>& buffer_dppd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForPPPD_hpp */
