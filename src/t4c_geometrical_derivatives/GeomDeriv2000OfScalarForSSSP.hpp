#ifndef GeomDeriv2000OfScalarForSSSP_hpp
#define GeomDeriv2000OfScalarForSSSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[SS|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_sssp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sssp: the primitive integrals buffer.
/// - Parameter buffer_dssp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_sssp_0(CSimdArray<double>& buffer_2000_sssp,
                     const CSimdArray<double>& buffer_sssp,
                     const CSimdArray<double>& buffer_dssp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForSSSP_hpp */
