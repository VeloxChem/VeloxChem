#ifndef NuclearPotentialGridFunc_hpp
#define NuclearPotentialGridFunc_hpp

#include "SubMatrix.hpp"

namespace npotfunc {

/// @brief Computes nuclear potential integrals on grid from pair of basis functions.
/// @param spher_buffer The spherical integrals buffer.
/// @param cart_buffer The Cartesian integrals buffer.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_igto The index of basis function on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
auto
compute(CSubMatrix&                      spher_buffer,
        CSubMatrix&                      cart_buffer,
        const CGtoBlock&                 bra_gto_block,
        const CGtoBlock&                 ket_gto_block,
        const int                        bra_igto,
        const int                        ket_igto,
        const bool                       bra_eq_ket) -> void
{
    const auto bra_angmom = bra_gto_block.angular_momentum();

    const auto ket_angmom = ket_gto_block.angular_momentum();

//    if ((bra_angmom == 0) && (ket_angmom == 0))
//    {
//        npotrec::comp_sum_nuclear_potential_ss(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);
//
//        return;
//    }
}

}  // namespace npotfunc

#endif /* NuclearPotentialGridFunc_hpp */
