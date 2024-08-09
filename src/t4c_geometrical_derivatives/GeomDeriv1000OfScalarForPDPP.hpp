#ifndef GeomDeriv1000OfScalarForPDPP_hpp
#define GeomDeriv1000OfScalarForPDPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[PD|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_pdpp: the integral geometrical derivatives buffer.
/// - Parameter buffer_sdpp: the primitive integrals buffer.
/// - Parameter buffer_ddpp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_pdpp_0(CSimdArray<double>& buffer_1000_pdpp,
                     const CSimdArray<double>& buffer_sdpp,
                     const CSimdArray<double>& buffer_ddpp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForPDPP_hpp */
