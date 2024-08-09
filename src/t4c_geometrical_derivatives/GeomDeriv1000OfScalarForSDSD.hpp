#ifndef GeomDeriv1000OfScalarForSDSD_hpp
#define GeomDeriv1000OfScalarForSDSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SD|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_sdsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdsd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_sdsd_0(CSimdArray<double>& buffer_1000_sdsd,
                     const CSimdArray<double>& buffer_pdsd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSDSD_hpp */
