#ifndef NuclearPotentialGridFunc_hpp
#define NuclearPotentialGridFunc_hpp

#include <vector>

#include "SubMatrix.hpp"
#include "NuclearPotentialGridRecSS.hpp"
#include "NuclearPotentialGridRecSP.hpp"
#include "NuclearPotentialGridRecPS.hpp"
#include "NuclearPotentialGridRecSD.hpp"
#include "NuclearPotentialGridRecDS.hpp"
#include "NuclearPotentialGridRecPP.hpp"
#include "NuclearPotentialGridRecSF.hpp"
#include "NuclearPotentialGridRecFS.hpp"
#include "NuclearPotentialGridRecPD.hpp"
#include "NuclearPotentialGridRecDP.hpp"
#include "NuclearPotentialGridRecPF.hpp"
#include "NuclearPotentialGridRecFP.hpp"
#include "NuclearPotentialGridRecDD.hpp"
#include "NuclearPotentialGridRecDF.hpp"
#include "NuclearPotentialGridRecFD.hpp"
#include "NuclearPotentialGridRecFF.hpp"

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
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        npotrec::comp_on_grid_nuclear_potential_sp(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        npotrec::comp_on_grid_nuclear_potential_ps(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        npotrec::comp_on_grid_nuclear_potential_sd(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        npotrec::comp_on_grid_nuclear_potential_ds(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        npotrec::comp_on_grid_nuclear_potential_pp(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        npotrec::comp_on_grid_nuclear_potential_sf(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        npotrec::comp_on_grid_nuclear_potential_fs(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        npotrec::comp_on_grid_nuclear_potential_pd(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        npotrec::comp_on_grid_nuclear_potential_dp(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        npotrec::comp_on_grid_nuclear_potential_pf(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        npotrec::comp_on_grid_nuclear_potential_fp(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        npotrec::comp_on_grid_nuclear_potential_dd(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        npotrec::comp_on_grid_nuclear_potential_df(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        npotrec::comp_on_grid_nuclear_potential_fd(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        npotrec::comp_on_grid_nuclear_potential_ff(spher_buffer, cart_buffer, gcoords_x, gcoords_y, gcoords_z, gweights, bra_gto_block, ket_gto_block, bra_igto, ket_igto);

        return;
    }
}

}  // namespace npotfunc

#endif /* NuclearPotentialGridFunc_hpp */
