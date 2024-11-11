#ifndef NuclearPotentialErfFunc_hpp
#define NuclearPotentialErfFunc_hpp

#include <array>

#include "GtoBlock.hpp"
#include "NuclearPotentialSumErfRecDD.hpp"
#include "NuclearPotentialSumErfRecDF.hpp"
#include "NuclearPotentialSumErfRecDG.hpp"
#include "NuclearPotentialSumErfRecDP.hpp"
#include "NuclearPotentialSumErfRecDS.hpp"
#include "NuclearPotentialSumErfRecFD.hpp"
#include "NuclearPotentialSumErfRecFF.hpp"
#include "NuclearPotentialSumErfRecFG.hpp"
#include "NuclearPotentialSumErfRecFP.hpp"
#include "NuclearPotentialSumErfRecFS.hpp"
#include "NuclearPotentialSumErfRecGD.hpp"
#include "NuclearPotentialSumErfRecGF.hpp"
#include "NuclearPotentialSumErfRecGG.hpp"
#include "NuclearPotentialSumErfRecGP.hpp"
#include "NuclearPotentialSumErfRecGS.hpp"
#include "NuclearPotentialSumErfRecPD.hpp"
#include "NuclearPotentialSumErfRecPF.hpp"
#include "NuclearPotentialSumErfRecPG.hpp"
#include "NuclearPotentialSumErfRecPP.hpp"
#include "NuclearPotentialSumErfRecPS.hpp"
#include "NuclearPotentialSumErfRecSD.hpp"
#include "NuclearPotentialSumErfRecSF.hpp"
#include "NuclearPotentialSumErfRecSG.hpp"
#include "NuclearPotentialSumErfRecSP.hpp"
#include "NuclearPotentialSumErfRecSS.hpp"

namespace npotfunc {

/// @brief Computes range separated nuclear potential integrals for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param omegas The vector of range-separation factors.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
compute(T&                               distributor,
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
        npotrec::comp_sum_erf_nuclear_potential_ss(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_sp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        npotrec::comp_sum_erf_nuclear_potential_ps(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_pp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_sd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        npotrec::comp_sum_erf_nuclear_potential_ds(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_sf(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        npotrec::comp_sum_erf_nuclear_potential_fs(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_pd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_dp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_sg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        npotrec::comp_sum_erf_nuclear_potential_gs(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_pf(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_fp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_dd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_pg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_gp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_df(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_fd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_dg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_gd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_ff(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_fg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_gf(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_gg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace npotfunc

#endif /* NuclearPotentialErfFunc_hpp */
