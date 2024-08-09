#ifndef GeomDeriv2000OfScalarForPPSD_hpp
#define GeomDeriv2000OfScalarForPPSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[PP|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_ppsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppsd: the primitive integrals buffer.
/// - Parameter buffer_fpsd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_ppsd_0(CSimdArray<double>& buffer_2000_ppsd,
                     const CSimdArray<double>& buffer_ppsd,
                     const CSimdArray<double>& buffer_fpsd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForPPSD_hpp */
