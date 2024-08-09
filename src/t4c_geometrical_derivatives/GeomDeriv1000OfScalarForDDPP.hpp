#ifndef GeomDeriv1000OfScalarForDDPP_hpp
#define GeomDeriv1000OfScalarForDDPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(1)/dA^(1)[DD|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_1000_ddpp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdpp: the primitive integrals buffer.
/// - Parameter buffer_fdpp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom1000_ddpp_0(CSimdArray<double>& buffer_1000_ddpp,
                     const CSimdArray<double>& buffer_pdpp,
                     const CSimdArray<double>& buffer_fdpp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv1000OfScalarForDDPP_hpp */
