#ifndef GeomDeriv1010OfScalarForPDPP_hpp
#define GeomDeriv1010OfScalarForPDPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[PD|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_pdpp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdsp: the primitive integrals buffer.
/// - Parameter buffer_sddp: the primitive integrals buffer.
/// - Parameter buffer_ddsp: the primitive integrals buffer.
/// - Parameter buffer_dddp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_pdpp_0(CSimdArray<double>& buffer_1010_pdpp,
                     const CSimdArray<double>& buffer_sdsp,
                     const CSimdArray<double>& buffer_sddp,
                     const CSimdArray<double>& buffer_ddsp,
                     const CSimdArray<double>& buffer_dddp,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForPDPP_hpp */
