#ifndef GeomDeriv1000OfScalarForPPSD_hpp
#define GeomDeriv1000OfScalarForPPSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[PP|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_ppsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_spsd: the primitive integrals buffer.
/// - Parameter buffer_dpsd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_ppsd_0(CSimdArray<double>& buffer_1000_ppsd,
                     const CSimdArray<double>& buffer_spsd,
                     const CSimdArray<double>& buffer_dpsd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForPPSD_hpp */
