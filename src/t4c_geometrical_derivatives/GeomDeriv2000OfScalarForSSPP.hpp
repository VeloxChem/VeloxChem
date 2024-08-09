#ifndef GeomDeriv2000OfScalarForSSPP_hpp
#define GeomDeriv2000OfScalarForSSPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SS|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_sspp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sspp: the primitive integrals buffer.
/// - Parameter buffer_dspp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_sspp_0(CSimdArray<double>& buffer_2000_sspp,
                     const CSimdArray<double>& buffer_sspp,
                     const CSimdArray<double>& buffer_dspp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSSPP_hpp */
