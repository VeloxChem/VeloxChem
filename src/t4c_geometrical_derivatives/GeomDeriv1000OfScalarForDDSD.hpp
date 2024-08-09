#ifndef GeomDeriv1000OfScalarForDDSD_hpp
#define GeomDeriv1000OfScalarForDDSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[DD|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_ddsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdsd: the primitive integrals buffer.
/// - Parameter buffer_fdsd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_ddsd_0(CSimdArray<double>& buffer_1000_ddsd,
                     const CSimdArray<double>& buffer_pdsd,
                     const CSimdArray<double>& buffer_fdsd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForDDSD_hpp */
