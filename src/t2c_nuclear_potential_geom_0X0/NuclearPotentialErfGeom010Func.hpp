#ifndef NuclearPotentialErfGeom010Func_hpp
#define NuclearPotentialErfGeom010Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "NuclearPotentialGeom010SumErfRecDD.hpp"
#include "NuclearPotentialGeom010SumErfRecDF.hpp"
#include "NuclearPotentialGeom010SumErfRecDG.hpp"
#include "NuclearPotentialGeom010SumErfRecDP.hpp"
#include "NuclearPotentialGeom010SumErfRecDS.hpp"
#include "NuclearPotentialGeom010SumErfRecFD.hpp"
#include "NuclearPotentialGeom010SumErfRecFF.hpp"
#include "NuclearPotentialGeom010SumErfRecFG.hpp"
#include "NuclearPotentialGeom010SumErfRecFP.hpp"
#include "NuclearPotentialGeom010SumErfRecFS.hpp"
#include "NuclearPotentialGeom010SumErfRecGD.hpp"
#include "NuclearPotentialGeom010SumErfRecGF.hpp"
#include "NuclearPotentialGeom010SumErfRecGG.hpp"
#include "NuclearPotentialGeom010SumErfRecGP.hpp"
#include "NuclearPotentialGeom010SumErfRecGS.hpp"
#include "NuclearPotentialGeom010SumErfRecPD.hpp"
#include "NuclearPotentialGeom010SumErfRecPF.hpp"
#include "NuclearPotentialGeom010SumErfRecPG.hpp"
#include "NuclearPotentialGeom010SumErfRecPP.hpp"
#include "NuclearPotentialGeom010SumErfRecPS.hpp"
#include "NuclearPotentialGeom010SumErfRecSD.hpp"
#include "NuclearPotentialGeom010SumErfRecSF.hpp"
#include "NuclearPotentialGeom010SumErfRecSG.hpp"
#include "NuclearPotentialGeom010SumErfRecSP.hpp"
#include "NuclearPotentialGeom010SumErfRecSS.hpp"

namespace npotfunc {

/// @brief Computes electric field integrals for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param omegas The vector of range-separation factors.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
compute_geom_010(T&                               distributor,
                 const std::vector<double>&       omegas,
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
        npotrec::comp_sum_erf_nuclear_potential_geom_010_ss(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_sp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_ps(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_pp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_sd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_ds(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_sf(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_fs(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_pd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_dp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_sg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_gs(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_pf(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_fp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_dd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_pg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_gp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_df(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_fd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_dg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_gd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_ff(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_fg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_gf(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_010_gg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace npotfunc

#endif /* NuclearPotentialErfGeom010Func_hpp */
