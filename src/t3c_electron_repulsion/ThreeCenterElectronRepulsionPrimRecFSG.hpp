#ifndef ThreeCenterElectronRepulsionPrimRecFSG_hpp
#define ThreeCenterElectronRepulsionPrimRecFSG_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [F|1/|r-r'||SG]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_fsg The index of integral in primitive integrals buffer.
/// @param idx_eri_0_psg The primitive integrals buffer.
/// @param idx_eri_1_psg The primitive integrals buffer.
/// @param idx_eri_1_dsf The primitive integrals buffer.
/// @param idx_eri_1_dsg The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_wa The vector of distances R(WA) = W - A.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_fsg(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_fsg,
                                 size_t idx_eri_0_psg,
                                 size_t idx_eri_1_psg,
                                 size_t idx_eri_1_dsf,
                                 size_t idx_eri_1_dsg,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecFSG_hpp */
