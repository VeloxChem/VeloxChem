#ifndef GeomDeriv1000OfScalarForDDSP_hpp
#define GeomDeriv1000OfScalarForDDSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[DD|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_ddsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdsp: the primitive integrals buffer.
/// - Parameter buffer_fdsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_ddsp_0(CSimdArray<double>& buffer_1000_ddsp,
                     const CSimdArray<double>& buffer_pdsp,
                     const CSimdArray<double>& buffer_fdsp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForDDSP_hpp */
