#ifndef GeomDeriv2000OfScalarForSDPD_hpp
#define GeomDeriv2000OfScalarForSDPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SD|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_sdpd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdpd: the primitive integrals buffer.
/// - Parameter buffer_ddpd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_sdpd_0(CSimdArray<double>& buffer_2000_sdpd,
                     const CSimdArray<double>& buffer_sdpd,
                     const CSimdArray<double>& buffer_ddpd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSDPD_hpp */
