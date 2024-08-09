#ifndef GeomDeriv1010OfScalarForSDSS_hpp
#define GeomDeriv1010OfScalarForSDSS_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[SD|G|SS]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_sdss: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdps: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_sdss_0(CSimdArray<double>& buffer_1010_sdss,
                     const CSimdArray<double>& buffer_pdps,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForSDSS_hpp */
