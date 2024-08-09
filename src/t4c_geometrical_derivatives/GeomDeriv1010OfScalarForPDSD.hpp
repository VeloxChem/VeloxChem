#ifndef GeomDeriv1010OfScalarForPDSD_hpp
#define GeomDeriv1010OfScalarForPDSD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)d^(1)/dC^(1)[PD|G|SD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1010_pdsd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdpd: the primitive integrals buffer.
/// - Parameter buffer_ddpd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
/// - Parameter c_exps: the vector of exponents on center C.
auto
comp_geom1010_pdsd_0(CSimdArray<double>& buffer_1010_pdsd,
                     const CSimdArray<double>& buffer_sdpd,
                     const CSimdArray<double>& buffer_ddpd,
                     const double a_exp,
                     const double* c_exps) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1010OfScalarForPDSD_hpp */
