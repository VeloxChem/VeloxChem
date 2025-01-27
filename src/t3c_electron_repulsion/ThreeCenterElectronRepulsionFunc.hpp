#ifndef ThreeCenterElectronRepulsionFunc_hpp
#define ThreeCenterElectronRepulsionFunc_hpp

#include <array>

#include "GtoBlock.hpp"
#include "GtoPairBlock.hpp"

#include "ThreeCenterElectronRepulsionRecSSS.hpp"
#include "ThreeCenterElectronRepulsionRecSSP.hpp"
#include "ThreeCenterElectronRepulsionRecSSD.hpp"
#include "ThreeCenterElectronRepulsionRecSPP.hpp"
#include "ThreeCenterElectronRepulsionRecSPD.hpp"
#include "ThreeCenterElectronRepulsionRecSDD.hpp"

#include "ThreeCenterElectronRepulsionRecPSS.hpp"
#include "ThreeCenterElectronRepulsionRecPSP.hpp"
#include "ThreeCenterElectronRepulsionRecPSD.hpp"
#include "ThreeCenterElectronRepulsionRecPPP.hpp"
#include "ThreeCenterElectronRepulsionRecPPD.hpp"
#include "ThreeCenterElectronRepulsionRecPDD.hpp"

#include "ThreeCenterElectronRepulsionRecDSS.hpp"
#include "ThreeCenterElectronRepulsionRecDSP.hpp"
#include "ThreeCenterElectronRepulsionRecDSD.hpp"
#include "ThreeCenterElectronRepulsionRecDPP.hpp"
#include "ThreeCenterElectronRepulsionRecDPD.hpp"
#include "ThreeCenterElectronRepulsionRecDDD.hpp"

namespace t3cerifunc {

/// @brief Computes electron repulsion integrals for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param aux_gto_block The basis functions block on auxilary side side.
/// @param gto_pair_block The basis function pairs block on ket side.
/// @param aux_indices The range [aux_first, aux_last) of basis functions on bra side.
template <class T>
auto
compute(T&                               distributor,
        const CGtoBlock&                 aux_gto_block,
        const CGtoPairBlock&             gto_pair_block,
        const std::pair<size_t, size_t>& aux_indices) -> void
{
    const auto aux_angmom = aux_gto_block.angular_momentum();

    const auto ket_angmoms = gto_pair_block.angular_momentums();

    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_sss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
         t3ceri::comp_electron_repulsion_ssp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
         t3ceri::comp_electron_repulsion_ssd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
         t3ceri::comp_electron_repulsion_spp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
         t3ceri::comp_electron_repulsion_spd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
         t3ceri::comp_electron_repulsion_sdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_pss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
         t3ceri::comp_electron_repulsion_psp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
         t3ceri::comp_electron_repulsion_psd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
         t3ceri::comp_electron_repulsion_ppp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
         t3ceri::comp_electron_repulsion_ppd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
         t3ceri::comp_electron_repulsion_pdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
   
    
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_dss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
         t3ceri::comp_electron_repulsion_dsp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
         t3ceri::comp_electron_repulsion_dsd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
         t3ceri::comp_electron_repulsion_dpp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
         t3ceri::comp_electron_repulsion_dpd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
         t3ceri::comp_electron_repulsion_ddd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
}

}  // namespace t2cerifunc


#endif /* ThreeCenterElectronRepulsionFunc_hpp */
