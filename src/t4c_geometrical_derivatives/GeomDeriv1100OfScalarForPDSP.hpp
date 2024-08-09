#ifndef GeomDeriv1100OfScalarForPDSP_hpp
#define GeomDeriv1100OfScalarForPDSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[PD|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_pdsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_spsp: the primitive integrals buffer.
/// - Parameter buffer_sfsp: the primitive integrals buffer.
/// - Parameter buffer_dpsp: the primitive integrals buffer.
/// - Parameter buffer_dfsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_pdsp_0(CSimdArray<double>& buffer_1100_pdsp,
                     const CSimdArray<double>& buffer_spsp,
                     const CSimdArray<double>& buffer_sfsp,
                     const CSimdArray<double>& buffer_dpsp,
                     const CSimdArray<double>& buffer_dfsp,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForPDSP_hpp */
