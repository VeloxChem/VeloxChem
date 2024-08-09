#ifndef GeomDeriv2000OfScalarForSSSS_hpp
#define GeomDeriv2000OfScalarForSSSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SS|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_ssss: the integral geometrical derivatives buffer.
/// - Parameter buffer_ssss: the primitive integrals buffer.
/// - Parameter buffer_dsss: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_ssss_0(CSimdArray<double>& buffer_2000_ssss,
                     const CSimdArray<double>& buffer_ssss,
                     const CSimdArray<double>& buffer_dsss,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSSSS_hpp */
