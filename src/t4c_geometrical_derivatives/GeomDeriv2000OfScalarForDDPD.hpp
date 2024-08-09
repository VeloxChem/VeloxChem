#ifndef GeomDeriv2000OfScalarForDDPD_hpp
#define GeomDeriv2000OfScalarForDDPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[DD|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_ddpd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdpd: the primitive integrals buffer.
/// - Parameter buffer_ddpd: the primitive integrals buffer.
/// - Parameter buffer_gdpd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_ddpd_0(CSimdArray<double>& buffer_2000_ddpd,
                     const CSimdArray<double>& buffer_sdpd,
                     const CSimdArray<double>& buffer_ddpd,
                     const CSimdArray<double>& buffer_gdpd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForDDPD_hpp */
