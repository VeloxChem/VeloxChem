#ifndef GeomDeriv2000OfScalarForPDSS_hpp
#define GeomDeriv2000OfScalarForPDSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[PD|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_pdss: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdss: the primitive integrals buffer.
/// - Parameter buffer_fdss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_pdss_0(CSimdArray<double>& buffer_2000_pdss,
                     const CSimdArray<double>& buffer_pdss,
                     const CSimdArray<double>& buffer_fdss,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForPDSS_hpp */
