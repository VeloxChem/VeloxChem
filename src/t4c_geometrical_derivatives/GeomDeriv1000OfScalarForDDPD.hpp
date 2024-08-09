#ifndef GeomDeriv1000OfScalarForDDPD_hpp
#define GeomDeriv1000OfScalarForDDPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[DD|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_ddpd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdpd: the primitive integrals buffer.
/// - Parameter buffer_fdpd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_ddpd_0(CSimdArray<double>& buffer_1000_ddpd,
                     const CSimdArray<double>& buffer_pdpd,
                     const CSimdArray<double>& buffer_fdpd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForDDPD_hpp */
