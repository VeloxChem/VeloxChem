#ifndef ThreeCenterElectronRepulsionPrimRecSSI_hpp
#define ThreeCenterElectronRepulsionPrimRecSSI_hpp

#include <cstddef>

#include "Point.hpp"
#include "SimdArray.hpp"

namespace t3ceri { // t3ceri namespace

/// Computes [S|1/|r-r'||SI]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_0_ssi The index of integral in primitive integrals buffer.
/// @param idx_eri_0_ssg The primitive integrals buffer.
/// @param idx_eri_1_ssg The primitive integrals buffer.
/// @param idx_eri_0_ssh The primitive integrals buffer.
/// @param idx_eri_1_ssh The primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_qd The vector of distances R(QD) = Q - D.
/// @param idx_wq The vector of distances R(WQ) = W - Q.
/// @param a_exp The exponent on center A.
auto
comp_prim_electron_repulsion_ssi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ssi,
                                 size_t idx_eri_0_ssg,
                                 size_t idx_eri_1_ssg,
                                 size_t idx_eri_0_ssh,
                                 size_t idx_eri_1_ssh,
                                 CSimdArray<double>& factors,
                                 const size_t idx_qd,
                                 const size_t idx_wq,
                                 const double a_exp) -> void;
} // t3ceri namespace

#endif /* ThreeCenterElectronRepulsionPrimRecSSI_hpp */
