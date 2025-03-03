#ifndef TwoCenterElectronRepulsionPrimRecSS_hpp
#define TwoCenterElectronRepulsionPrimRecSS_hpp

#include "SimdArray.hpp"

namespace t2ceri { // t2ceri namespace

/// @brief Computes primitive [S|1/|r-r'||S]  integrals for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_eri_ss The index of integral in primitive integrals buffer.

auto
comp_prim_electron_repulsion_ss(CSimdArray<double>& pbuffer,
                                const size_t idx_eri_ss,
                                const CSimdArray<double>& bf_data,
                                const size_t              idx_vals,
                                CSimdArray<double>&       factors,
                                const double              a_exp,
                                const double              a_norm) -> void;
} // t2ceri namespace

#endif /* TwoCenterElectronRepulsionPrimRecSS_hpp */
