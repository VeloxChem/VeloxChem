#ifndef GeomDeriv2000OfScalarForDDSP_hpp
#define GeomDeriv2000OfScalarForDDSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[DD|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_ddsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdsp: the primitive integrals buffer.
/// - Parameter buffer_ddsp: the primitive integrals buffer.
/// - Parameter buffer_gdsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_ddsp_0(CSimdArray<double>& buffer_2000_ddsp,
                     const CSimdArray<double>& buffer_sdsp,
                     const CSimdArray<double>& buffer_ddsp,
                     const CSimdArray<double>& buffer_gdsp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForDDSP_hpp */
