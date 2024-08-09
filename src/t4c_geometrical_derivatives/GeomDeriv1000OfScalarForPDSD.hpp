#ifndef GeomDeriv1000OfScalarForPDSD_hpp
#define GeomDeriv1000OfScalarForPDSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[PD|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_pdsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdsd: the primitive integrals buffer.
/// - Parameter buffer_ddsd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_pdsd_0(CSimdArray<double>& buffer_1000_pdsd,
                     const CSimdArray<double>& buffer_sdsd,
                     const CSimdArray<double>& buffer_ddsd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForPDSD_hpp */
