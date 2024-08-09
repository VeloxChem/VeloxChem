#ifndef GeomDeriv1010OfScalarForDDDD_hpp
#define GeomDeriv1010OfScalarForDDDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[DD|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_dddd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdpd: the primitive integrals buffer.
/// - Parameter buffer_pdfd: the primitive integrals buffer.
/// - Parameter buffer_fdpd: the primitive integrals buffer.
/// - Parameter buffer_fdfd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_dddd_0(CSimdArray<double>& buffer_1010_dddd,
                     const CSimdArray<double>& buffer_pdpd,
                     const CSimdArray<double>& buffer_pdfd,
                     const CSimdArray<double>& buffer_fdpd,
                     const CSimdArray<double>& buffer_fdfd,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForDDDD_hpp */
