#ifndef GeomDeriv1000OfScalarForSSPP_hpp
#define GeomDeriv1000OfScalarForSSPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SS|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_sspp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pspp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_sspp_0(CSimdArray<double>& buffer_1000_sspp,
                     const CSimdArray<double>& buffer_pspp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSSPP_hpp */
