#ifndef TwoCenterElectronRepulsionPrimRecPG
#define TwoCenterElectronRepulsionPrimRecPG

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [P|1/|r-r'||G]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_pg The index of integral in primitive integrals buffer.
/// @param idx_eri_1_sf The index of integral in primitive integrals buffer.
/// @param idx_eri_1_sg The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_rpa The vector of distances R(PA) = P - A.
/// @param a_exp The primitive basis function exponent on center A.
auto
comp_prim_electron_repulsion_pg(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_pg,
                                const size_t idx_eri_1_sf,
                                const size_t idx_eri_1_sg,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecPG */
