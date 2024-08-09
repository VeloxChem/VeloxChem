#ifndef GeomDeriv1000OfScalarForSSDD_hpp
#define GeomDeriv1000OfScalarForSSDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[SS|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_ssdd: the integral geometrical derivatives buffer.
/// - Parameter buffer_psdd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_ssdd_0(CSimdArray<double>& buffer_1000_ssdd,
                     const CSimdArray<double>& buffer_psdd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForSSDD_hpp */
