#ifndef ThreeCenterElectronRepulsionGeom010Func_hpp
#define ThreeCenterElectronRepulsionGeom010Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "GtoPairBlock.hpp"

#include "ThreeCenterElectronRepulsionGeom010RecSSS.hpp"

namespace t3cerifunc {

/// @brief Computes electron repulsion integrals for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param aux_gto_block The basis functions block on auxilary side side.
/// @param gto_pair_block The basis function pairs block on ket side.
/// @param aux_indices The range [aux_first, aux_last) of basis functions on bra side.
template <class T>
auto
compute_geom_010(T&                               distributor,
                 const CGtoBlock&                 aux_gto_block,
                 const CGtoPairBlock&             gto_pair_block,
                 const std::pair<size_t, size_t>& aux_indices) -> void
{
    const auto aux_angmom = aux_gto_block.angular_momentum();

    const auto ket_angmoms = gto_pair_block.angular_momentums();

    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_sss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
}

}  // namespace t2cerifunc


#endif /* ThreeCenterElectronRepulsionGeom010Func_hpp */
