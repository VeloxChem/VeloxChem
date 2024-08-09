#ifndef GeomDeriv1000OfScalarForPDSS_hpp
#define GeomDeriv1000OfScalarForPDSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[PD|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_pdss: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdss: the primitive integrals buffer.
/// - Parameter buffer_ddss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_pdss_0(CSimdArray<double>& buffer_1000_pdss,
                     const CSimdArray<double>& buffer_sdss,
                     const CSimdArray<double>& buffer_ddss,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForPDSS_hpp */
