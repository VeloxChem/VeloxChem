#ifndef GeomDeriv1000OfScalarForSSSD_hpp
#define GeomDeriv1000OfScalarForSSSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SS|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_sssd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pssd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_sssd_0(CSimdArray<double>& buffer_1000_sssd,
                     const CSimdArray<double>& buffer_pssd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSSSD_hpp */
