#ifndef GeomDeriv1000OfScalarForDDSS_hpp
#define GeomDeriv1000OfScalarForDDSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[DD|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_ddss: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdss: the primitive integrals buffer.
/// - Parameter buffer_fdss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_ddss_0(CSimdArray<double>& buffer_1000_ddss,
                     const CSimdArray<double>& buffer_pdss,
                     const CSimdArray<double>& buffer_fdss,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForDDSS_hpp */
