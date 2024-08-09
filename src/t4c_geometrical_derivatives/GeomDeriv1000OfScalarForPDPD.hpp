#ifndef GeomDeriv1000OfScalarForPDPD_hpp
#define GeomDeriv1000OfScalarForPDPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[PD|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_pdpd: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdpd: the primitive integrals buffer.
/// - Parameter buffer_ddpd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_pdpd_0(CSimdArray<double>& buffer_1000_pdpd,
                     const CSimdArray<double>& buffer_sdpd,
                     const CSimdArray<double>& buffer_ddpd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForPDPD_hpp */
