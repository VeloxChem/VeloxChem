#ifndef GeomDeriv1010OfScalarForDDSP_hpp
#define GeomDeriv1010OfScalarForDDSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[DD|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_ddsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdpp: the primitive integrals buffer.
/// - Parameter buffer_fdpp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_ddsp_0(CSimdArray<double>& buffer_1010_ddsp,
                     const CSimdArray<double>& buffer_pdpp,
                     const CSimdArray<double>& buffer_fdpp,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForDDSP_hpp */
