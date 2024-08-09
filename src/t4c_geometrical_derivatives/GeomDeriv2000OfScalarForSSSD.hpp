#ifndef GeomDeriv2000OfScalarForSSSD_hpp
#define GeomDeriv2000OfScalarForSSSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SS|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_sssd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sssd: the primitive integrals buffer.
/// - Parameter buffer_dssd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_sssd_0(CSimdArray<double>& buffer_2000_sssd,
                     const CSimdArray<double>& buffer_sssd,
                     const CSimdArray<double>& buffer_dssd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSSSD_hpp */
