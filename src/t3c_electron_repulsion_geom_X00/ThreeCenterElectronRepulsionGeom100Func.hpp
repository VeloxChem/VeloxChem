#ifndef ThreeCenterElectronRepulsionGeom100Func_hpp
#define ThreeCenterElectronRepulsionGeom100Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "GtoPairBlock.hpp"

#include "ThreeCenterElectronRepulsionGeom100RecSSS.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSSP.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSSD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSPP.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSSF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSPD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSSG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSPF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSDD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSPG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSDF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSDG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSFF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSFG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSGG.hpp"

namespace t3cerifunc {

/// @brief Computes electron repulsion integrals for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param aux_gto_block The basis functions block on auxilary side side.
/// @param gto_pair_block The basis function pairs block on ket side.
/// @param aux_indices The range [aux_first, aux_last) of basis functions on bra side.
template <class T>
auto
compute_geom_100(T&                               distributor,
                 const CGtoBlock&                 aux_gto_block,
                 const CGtoPairBlock&             gto_pair_block,
                 const std::pair<size_t, size_t>& aux_indices) -> void
{
    const auto aux_angmom = aux_gto_block.angular_momentum();

    const auto ket_angmoms = gto_pair_block.angular_momentums();

    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_geom100_sss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        t3ceri::comp_electron_repulsion_geom100_ssp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_ssd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        t3ceri::comp_electron_repulsion_geom100_spp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_ssf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_spd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_ssg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_spf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_sdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_spg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_sdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_sdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_sff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_sfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_sgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
}

}  // namespace t2cerifunc

#endif /* ThreeCenterElectronRepulsionGeom100Func_hpp */
