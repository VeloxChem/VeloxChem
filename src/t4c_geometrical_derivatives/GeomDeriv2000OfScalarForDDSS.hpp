#ifndef GeomDeriv2000OfScalarForDDSS_hpp
#define GeomDeriv2000OfScalarForDDSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[DD|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_ddss: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdss: the primitive integrals buffer.
/// - Parameter buffer_ddss: the primitive integrals buffer.
/// - Parameter buffer_gdss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_ddss_0(CSimdArray<double>& buffer_2000_ddss,
                     const CSimdArray<double>& buffer_sdss,
                     const CSimdArray<double>& buffer_ddss,
                     const CSimdArray<double>& buffer_gdss,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForDDSS_hpp */
