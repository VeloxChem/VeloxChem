#ifndef NuclearPotentialGeom100Func_hpp
#define NuclearPotentialGeom100Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "NuclearPotentialGeom100SumRecDD.hpp"
#include "NuclearPotentialGeom100SumRecDF.hpp"
#include "NuclearPotentialGeom100SumRecDG.hpp"
#include "NuclearPotentialGeom100SumRecDP.hpp"
#include "NuclearPotentialGeom100SumRecDS.hpp"
#include "NuclearPotentialGeom100SumRecFD.hpp"
#include "NuclearPotentialGeom100SumRecFF.hpp"
#include "NuclearPotentialGeom100SumRecFG.hpp"
#include "NuclearPotentialGeom100SumRecFP.hpp"
#include "NuclearPotentialGeom100SumRecFS.hpp"
#include "NuclearPotentialGeom100SumRecGD.hpp"
#include "NuclearPotentialGeom100SumRecGF.hpp"
#include "NuclearPotentialGeom100SumRecGG.hpp"
#include "NuclearPotentialGeom100SumRecGP.hpp"
#include "NuclearPotentialGeom100SumRecGS.hpp"
#include "NuclearPotentialGeom100SumRecPD.hpp"
#include "NuclearPotentialGeom100SumRecPF.hpp"
#include "NuclearPotentialGeom100SumRecPG.hpp"
#include "NuclearPotentialGeom100SumRecPP.hpp"
#include "NuclearPotentialGeom100SumRecPS.hpp"
#include "NuclearPotentialGeom100SumRecSD.hpp"
#include "NuclearPotentialGeom100SumRecSF.hpp"
#include "NuclearPotentialGeom100SumRecSG.hpp"
#include "NuclearPotentialGeom100SumRecSP.hpp"
#include "NuclearPotentialGeom100SumRecSS.hpp"

namespace npotfunc {

/// @brief Computes nuclear potential integral derivative for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
compute_geom_100(T&                               distributor,
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
        npotrec::comp_sum_nuclear_potential_geom_10_ss(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_sp(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_ps(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_pp(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_sd(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_ds(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_sf(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_fs(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_pd(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_dp(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_sg(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_gs(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_pf(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_fp(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_dd(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_pg(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_gp(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_df(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_fd(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_dg(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_gd(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_ff(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_fg(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_gf(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        npotrec::comp_sum_nuclear_potential_geom_10_gg(
            distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace npotfunc

#endif /* NuclearPotentialGeom100Func_hpp */
