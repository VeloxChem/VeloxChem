#ifndef GeomDeriv1010OfScalarForDDPD_hpp
#define GeomDeriv1010OfScalarForDDPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[DD|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_ddpd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdsd: the primitive integrals buffer.
/// - Parameter buffer_pddd: the primitive integrals buffer.
/// - Parameter buffer_fdsd: the primitive integrals buffer.
/// - Parameter buffer_fddd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_ddpd_0(CSimdArray<double>& buffer_1010_ddpd,
                     const CSimdArray<double>& buffer_pdsd,
                     const CSimdArray<double>& buffer_pddd,
                     const CSimdArray<double>& buffer_fdsd,
                     const CSimdArray<double>& buffer_fddd,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForDDPD_hpp */