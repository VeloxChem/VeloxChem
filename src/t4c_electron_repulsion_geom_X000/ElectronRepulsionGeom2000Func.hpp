#ifndef ElectronRepulsionGeom2000Func_hpp
#define ElectronRepulsionGeom2000Func_hpp

#include <cstddef>
#include <iostream>
#include <utility>

#include "ElectronRepulsionGeom2000RecSSSS.hpp"

#include "GtoPairBlock.hpp"

namespace erifunc {  // erifunc namespace

/// Computes vector of integrals of requested four center integral.
/// @param distributor  The pointer to distributor of integrals.
/// @param bra_gto_pair_block The basis function pairs block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
compute_geom_2000(T&                               distributor,
                  const CGtoPairBlock&             bra_gto_pair_block,
                  const CGtoPairBlock&             ket_gto_pair_block,
                  const std::pair<size_t, size_t>& bra_indices,
                  const std::pair<size_t, size_t>& ket_indices) -> void
{
    const auto bra_angmoms = bra_gto_pair_block.angular_momentums();

    const auto ket_angmoms = ket_gto_pair_block.angular_momentums();
    
    // leading [SS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom2000_ssss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    std::cout << " *** Integral not found in call tree :" << bra_angmoms.first << " , " << bra_angmoms.second << " , " << ket_angmoms.first << " , " << ket_angmoms.second << std::endl;
}

}  // namespace erifunc

#endif /* ElectronRepulsionGeom2000Func_hpp */
