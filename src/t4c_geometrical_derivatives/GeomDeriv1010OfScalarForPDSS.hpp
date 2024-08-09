#ifndef GeomDeriv1010OfScalarForPDSS_hpp
#define GeomDeriv1010OfScalarForPDSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[PD|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_pdss: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdps: the primitive integrals buffer.
/// - Parameter buffer_ddps: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_pdss_0(CSimdArray<double>& buffer_1010_pdss,
                     const CSimdArray<double>& buffer_sdps,
                     const CSimdArray<double>& buffer_ddps,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForPDSS_hpp */
