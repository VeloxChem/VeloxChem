#ifndef GeomDeriv1000OfScalarForSSSS_hpp
#define GeomDeriv1000OfScalarForSSSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SS|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_ssss: the integral geometrical derivatives buffer.
/// - Parameter buffer_psss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_ssss_0(CSimdArray<double>& buffer_1000_ssss,
                     const CSimdArray<double>& buffer_psss,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSSSS_hpp */
