#ifndef GeomDeriv1010OfScalarForPDPD_hpp
#define GeomDeriv1010OfScalarForPDPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[PD|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_pdpd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdsd: the primitive integrals buffer.
/// - Parameter buffer_sddd: the primitive integrals buffer.
/// - Parameter buffer_ddsd: the primitive integrals buffer.
/// - Parameter buffer_dddd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_pdpd_0(CSimdArray<double>& buffer_1010_pdpd,
                     const CSimdArray<double>& buffer_sdsd,
                     const CSimdArray<double>& buffer_sddd,
                     const CSimdArray<double>& buffer_ddsd,
                     const CSimdArray<double>& buffer_dddd,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForPDPD_hpp */
