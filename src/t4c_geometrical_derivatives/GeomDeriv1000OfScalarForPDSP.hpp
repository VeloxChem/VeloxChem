#ifndef GeomDeriv1000OfScalarForPDSP_hpp
#define GeomDeriv1000OfScalarForPDSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[PD|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_pdsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdsp: the primitive integrals buffer.
/// - Parameter buffer_ddsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_pdsp_0(CSimdArray<double>& buffer_1000_pdsp,
                     const CSimdArray<double>& buffer_sdsp,
                     const CSimdArray<double>& buffer_ddsp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForPDSP_hpp */
