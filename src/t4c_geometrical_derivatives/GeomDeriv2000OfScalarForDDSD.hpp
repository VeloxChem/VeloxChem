#ifndef GeomDeriv2000OfScalarForDDSD_hpp
#define GeomDeriv2000OfScalarForDDSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[DD|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_ddsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdsd: the primitive integrals buffer.
/// - Parameter buffer_ddsd: the primitive integrals buffer.
/// - Parameter buffer_gdsd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_ddsd_0(CSimdArray<double>& buffer_2000_ddsd,
                     const CSimdArray<double>& buffer_sdsd,
                     const CSimdArray<double>& buffer_ddsd,
                     const CSimdArray<double>& buffer_gdsd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForDDSD_hpp */
