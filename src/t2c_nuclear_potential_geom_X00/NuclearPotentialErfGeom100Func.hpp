//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef NuclearPotentialErfGeom100Func_hpp
#define NuclearPotentialErfGeom100Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "NuclearPotentialGeom100SumErfRecDD.hpp"
#include "NuclearPotentialGeom100SumErfRecDF.hpp"
#include "NuclearPotentialGeom100SumErfRecDG.hpp"
#include "NuclearPotentialGeom100SumErfRecDP.hpp"
#include "NuclearPotentialGeom100SumErfRecDS.hpp"
#include "NuclearPotentialGeom100SumErfRecFD.hpp"
#include "NuclearPotentialGeom100SumErfRecFF.hpp"
#include "NuclearPotentialGeom100SumErfRecFG.hpp"
#include "NuclearPotentialGeom100SumErfRecFP.hpp"
#include "NuclearPotentialGeom100SumErfRecFS.hpp"
#include "NuclearPotentialGeom100SumErfRecGD.hpp"
#include "NuclearPotentialGeom100SumErfRecGF.hpp"
#include "NuclearPotentialGeom100SumErfRecGG.hpp"
#include "NuclearPotentialGeom100SumErfRecGP.hpp"
#include "NuclearPotentialGeom100SumErfRecGS.hpp"
#include "NuclearPotentialGeom100SumErfRecPD.hpp"
#include "NuclearPotentialGeom100SumErfRecPF.hpp"
#include "NuclearPotentialGeom100SumErfRecPG.hpp"
#include "NuclearPotentialGeom100SumErfRecPP.hpp"
#include "NuclearPotentialGeom100SumErfRecPS.hpp"
#include "NuclearPotentialGeom100SumErfRecSD.hpp"
#include "NuclearPotentialGeom100SumErfRecSF.hpp"
#include "NuclearPotentialGeom100SumErfRecSG.hpp"
#include "NuclearPotentialGeom100SumErfRecSP.hpp"
#include "NuclearPotentialGeom100SumErfRecSS.hpp"

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
        npotrec::comp_sum_erf_nuclear_potential_geom_10_ss(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_sp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_ps(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_pp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_sd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_ds(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_sf(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_fs(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_pd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_dp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_sg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_gs(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_pf(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_fp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_dd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_pg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_gp(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_df(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_fd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_dg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_gd(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_ff(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_fg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_gf(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        npotrec::comp_sum_erf_nuclear_potential_geom_10_gg(distributor, omegas, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace npotfunc

#endif /* NuclearPotentialErfGeom100Func_hpp */
