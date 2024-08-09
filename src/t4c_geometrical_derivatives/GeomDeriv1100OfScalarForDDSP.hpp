#ifndef GeomDeriv1100OfScalarForDDSP_hpp
#define GeomDeriv1100OfScalarForDDSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dB^(1)[DD|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1100_ddsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_ppsp: the primitive integrals buffer.
/// - Parameter buffer_pfsp: the primitive integrals buffer.
/// - Parameter buffer_fpsp: the primitive integrals buffer.
/// - Parameter buffer_ffsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter b_exp: the exponent on center B.
auto
comp_geom1100_ddsp_0(CSimdArray<double>& buffer_1100_ddsp,
                     const CSimdArray<double>& buffer_ppsp,
                     const CSimdArray<double>& buffer_pfsp,
                     const CSimdArray<double>& buffer_fpsp,
                     const CSimdArray<double>& buffer_ffsp,
                     const double a_exp,
                     const double b_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1100OfScalarForDDSP_hpp */
