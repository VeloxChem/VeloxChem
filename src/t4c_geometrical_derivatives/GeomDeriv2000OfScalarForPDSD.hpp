#ifndef GeomDeriv2000OfScalarForPDSD_hpp
#define GeomDeriv2000OfScalarForPDSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[PD|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_pdsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdsd: the primitive integrals buffer.
/// - Parameter buffer_fdsd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_pdsd_0(CSimdArray<double>& buffer_2000_pdsd,
                     const CSimdArray<double>& buffer_pdsd,
                     const CSimdArray<double>& buffer_fdsd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForPDSD_hpp */
