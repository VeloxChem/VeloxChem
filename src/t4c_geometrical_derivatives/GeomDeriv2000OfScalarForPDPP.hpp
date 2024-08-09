#ifndef GeomDeriv2000OfScalarForPDPP_hpp
#define GeomDeriv2000OfScalarForPDPP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[PD|G|PP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_pdpp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdpp: the primitive integrals buffer.
/// - Parameter buffer_fdpp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_pdpp_0(CSimdArray<double>& buffer_2000_pdpp,
                     const CSimdArray<double>& buffer_pdpp,
                     const CSimdArray<double>& buffer_fdpp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForPDPP_hpp */
