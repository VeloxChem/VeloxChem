#ifndef ElectronRepulsionFunc_hpp
#define ElectronRepulsionFunc_hpp

#include <cstddef>
#include <utility>

#include "GtoPairBlock.hpp"

#include "ElectronRepulsionRecSSSS.hpp"
#include "ElectronRepulsionRecSSSP.hpp"
#include "ElectronRepulsionRecSSPP.hpp"
#include "ElectronRepulsionRecSPSS.hpp"
#include "ElectronRepulsionRecSPSP.hpp"
#include "ElectronRepulsionRecSPPP.hpp"
#include "ElectronRepulsionRecPPSS.hpp"
#include "ElectronRepulsionRecPPSP.hpp"
#include "ElectronRepulsionRecPPPP.hpp"

#include <iostream>

namespace erifunc {  // erifunc namespace

/// Computes vector of integrals of requested four center integral.
/// @param distributor  The pointer to distributor of integrals.
/// @param bra_gto_pair_block The basis function pairs block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
/// @param bra_eq_ket True if basis function pairs blocks on bra and ket are the same, False otherwise.
template <class T>
inline auto
compute(      T&                  distributor,
        const CGtoPairBlock& bra_gto_pair_block,
        const CGtoPairBlock& ket_gto_pair_block,
        const std::pair<size_t, size_t>& bra_indices,
        const std::pair<size_t, size_t>& ket_indices,
        const bool bra_eq_ket) -> void
{
    const auto bra_angmoms = bra_gto_pair_block.angular_momentums();
    
    const auto ket_angmoms = ket_gto_pair_block.angular_momentums();

    if ((bra_angmoms == std::pair<int, int>({0, 0})) &&
        (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_ssss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
            return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) &&
        (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_sssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
            return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 0})) &&
        (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_sspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
            return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) &&
        (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_spss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
            return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) &&
        (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_spsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
            return;
    }

    if ((bra_angmoms == std::pair<int, int>({0, 1})) &&
        (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_sppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
            return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) &&
        (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_ppss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices, bra_eq_ket);
            return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) &&
        (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_ppsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices,  bra_eq_ket);
            return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) &&
        (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_pppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices,  bra_eq_ket);
            return;
    }

    std::cout << " *** Integral not found in call tree :" << bra_angmoms.first << " , " << bra_angmoms.second  << " , " << ket_angmoms.first << " , " << ket_angmoms.second << std::endl;
}

}  // namespace erifunc

#endif /* ElectronRepulsionFunc_hpp */
