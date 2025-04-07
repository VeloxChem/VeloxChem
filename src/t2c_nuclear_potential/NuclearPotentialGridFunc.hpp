#ifndef NuclearPotentialGridFunc_hpp
#define NuclearPotentialGridFunc_hpp

#include <vector>

#include "SubMatrix.hpp"
#include "NuclearPotentialGridRecSS.hpp"

namespace npotfunc {

/// @brief Computes nuclear potential integrals on grid from pair of basis functions.
/// @param spher_buffer The spherical integrals buffer.
/// @param cart_buffer The Cartesian integrals buffer.
/// @param gcoords_x The Cartesian X coordinates of grid points.
/// @param gcoords_y The Cartesian Y coordinates of grid points.
/// @param gcoords_z The Cartesian Z coordinates of grid points.
/// @param gweights The weight of grid points.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_igto The index of basis function on ket side.
auto
compute(CSubMatrix&                spher_buffer,
        CSubMatrix&                cart_buffer,
        const std::vector<double>& gcoords_x,
        const std::vector<double>& gcoords_y,
        const std::vector<double>& gcoords_z,
        const std::vector<double>& gweights,
        const CGtoBlock&           bra_gto_block,
        const CGtoBlock&           ket_gto_block,
        const int                  bra_igto,
        const int                  ket_igto) -> void
{
    const auto bra_angmom = bra_gto_block.angular_momentum();

    const auto ket_angmom = ket_gto_block.angular_momentum();

    if ((bra_angmom == 0) && (ket_angmom == 0))
    {
        npotrec::comp_on_grid_nuclear_potential_ss(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
}

}  // namespace npotfunc

#endif /* NuclearPotentialGridFunc_hpp */
