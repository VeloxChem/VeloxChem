#ifndef GeomDeriv1000OfScalarForSSPD_hpp
#define GeomDeriv1000OfScalarForSSPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SS|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_sspd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pspd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_sspd_0(CSimdArray<double>& buffer_1000_sspd,
                     const CSimdArray<double>& buffer_pspd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSSPD_hpp */
