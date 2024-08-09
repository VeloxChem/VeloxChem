#ifndef GeomDeriv2000OfScalarForSSPD_hpp
#define GeomDeriv2000OfScalarForSSPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SS|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_sspd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sspd: the primitive integrals buffer.
/// - Parameter buffer_dspd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_sspd_0(CSimdArray<double>& buffer_2000_sspd,
                     const CSimdArray<double>& buffer_sspd,
                     const CSimdArray<double>& buffer_dspd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSSPD_hpp */
