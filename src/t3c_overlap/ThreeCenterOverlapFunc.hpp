#ifndef ThreeCenterOverlapFunc_hpp
#define ThreeCenterOverlapFunc_hpp

#include <array>

#include "GtoBlock.hpp"
#include "ThreeCenterOverlapRecSS.hpp"

namespace t3ovlfunc {

/// @brief Computes overlap integrals for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
compute(T&                               distributor,
        const CGtoBlock&                 bra_gto_block,
        const CGtoBlock&                 ket_gto_block,
        const std::pair<size_t, size_t>& bra_indices,
        const std::pair<size_t, size_t>& ket_indices,
        const bool                       bra_eq_ket) -> void
{
    const auto bra_angmom = bra_gto_block.angular_momentum();

    const auto ket_angmom = ket_gto_block.angular_momentum();

    if ((bra_angmom == 0) && (ket_angmom == 0))
    {
        t3ovlrec::comp_overlap_ss(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t3ovlfunc

#endif /* ThreeCenterOverlapFunc_hpp */
