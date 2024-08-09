#ifndef GeomDeriv1010OfScalarForPDDD_hpp
#define GeomDeriv1010OfScalarForPDDD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[PD|G|DD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_pddd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdpd: the primitive integrals buffer.
/// - Parameter buffer_sdfd: the primitive integrals buffer.
/// - Parameter buffer_ddpd: the primitive integrals buffer.
/// - Parameter buffer_ddfd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_pddd_0(CSimdArray<double>& buffer_1010_pddd,
                     const CSimdArray<double>& buffer_sdpd,
                     const CSimdArray<double>& buffer_sdfd,
                     const CSimdArray<double>& buffer_ddpd,
                     const CSimdArray<double>& buffer_ddfd,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForPDDD_hpp */
