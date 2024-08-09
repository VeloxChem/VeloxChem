#ifndef GeomDeriv2000OfScalarForPDPD_hpp
#define GeomDeriv2000OfScalarForPDPD_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[PD|G|PD]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_pdpd: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdpd: the primitive integrals buffer.
/// - Parameter buffer_fdpd: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_pdpd_0(CSimdArray<double>& buffer_2000_pdpd,
                     const CSimdArray<double>& buffer_pdpd,
                     const CSimdArray<double>& buffer_fdpd,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForPDPD_hpp */
