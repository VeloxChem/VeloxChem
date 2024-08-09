#ifndef GeomDeriv2000OfScalarForSSDD_hpp
#define GeomDeriv2000OfScalarForSSDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SS|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_ssdd: the integral geometrical derivatives buffer.
/// - Parameter buffer_ssdd: the primitive integrals buffer.
/// - Parameter buffer_dsdd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_ssdd_0(CSimdArray<double>& buffer_2000_ssdd,
                     const CSimdArray<double>& buffer_ssdd,
                     const CSimdArray<double>& buffer_dsdd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSSDD_hpp */
