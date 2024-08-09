#ifndef GeomDeriv2000OfScalarForPDSP_hpp
#define GeomDeriv2000OfScalarForPDSP_hpp

#include "SimdArray.hpp"

namespace t4c_geom { // t4c_geom namespace

/// Computes d^(2)/dA^(2)[PD|G|SP]  integrals for arbitrary scalar operator G.
/// - Parameter buffer_2000_pdsp: the integral geometrical derivatives buffer.
/// - Parameter buffer_pdsp: the primitive integrals buffer.
/// - Parameter buffer_fdsp: the primitive integrals buffer.
/// - Parameter a_exp: the exponent on center A.
auto
comp_geom2000_pdsp_0(CSimdArray<double>& buffer_2000_pdsp,
                     const CSimdArray<double>& buffer_pdsp,
                     const CSimdArray<double>& buffer_fdsp,
                     const double a_exp) -> void;

} // t4c_geom namespace

#endif /* GeomDeriv2000OfScalarForPDSP_hpp */
