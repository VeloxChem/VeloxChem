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

#ifndef LocalCorePotentialGeom101Func_hpp
#define LocalCorePotentialGeom101Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "LocalCorePotentialGeom101DD.hpp"
#include "LocalCorePotentialGeom101DF.hpp"
#include "LocalCorePotentialGeom101DG.hpp"
#include "LocalCorePotentialGeom101DP.hpp"
#include "LocalCorePotentialGeom101DS.hpp"
#include "LocalCorePotentialGeom101FD.hpp"
#include "LocalCorePotentialGeom101FF.hpp"
#include "LocalCorePotentialGeom101FG.hpp"
#include "LocalCorePotentialGeom101FP.hpp"
#include "LocalCorePotentialGeom101FS.hpp"
#include "LocalCorePotentialGeom101GD.hpp"
#include "LocalCorePotentialGeom101GF.hpp"
#include "LocalCorePotentialGeom101GG.hpp"
#include "LocalCorePotentialGeom101GP.hpp"
#include "LocalCorePotentialGeom101GS.hpp"
#include "LocalCorePotentialGeom101PD.hpp"
#include "LocalCorePotentialGeom101PF.hpp"
#include "LocalCorePotentialGeom101PG.hpp"
#include "LocalCorePotentialGeom101PP.hpp"
#include "LocalCorePotentialGeom101PS.hpp"
#include "LocalCorePotentialGeom101SD.hpp"
#include "LocalCorePotentialGeom101SF.hpp"
#include "LocalCorePotentialGeom101SG.hpp"
#include "LocalCorePotentialGeom101SP.hpp"
#include "LocalCorePotentialGeom101SS.hpp"

namespace t2lecp {

/// @brief Computes local ECP integrals for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param ecp_potential The local ECP potential.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
compute_geom_101(T&                               distributor,
                 const CGtoBlock&                 bra_gto_block,
                 const CGtoBlock&                 ket_gto_block,
                 const CBaseCorePotential&        ecp_potential,
                 const std::pair<size_t, size_t>& bra_indices,
                 const std::pair<size_t, size_t>& ket_indices,
                 const bool                       bra_eq_ket) -> void
{
    const auto bra_angmom = bra_gto_block.angular_momentum();

    const auto ket_angmom = ket_gto_block.angular_momentum();

    if ((bra_angmom == 0) && (ket_angmom == 0))
    {
        t2lecp::comp_local_core_potential_geom_101_ss(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_101_sp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        t2lecp::comp_local_core_potential_geom_101_ps(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_101_pp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_101_sd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        t2lecp::comp_local_core_potential_geom_101_ds(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_101_sf(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        t2lecp::comp_local_core_potential_geom_101_fs(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_101_pd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_101_dp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_101_sg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        t2lecp::comp_local_core_potential_geom_101_gs(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_101_pf(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_101_fp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_101_dd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_101_pg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_101_gp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_101_df(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_101_fd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_101_dg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_101_gd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_101_ff(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_101_fg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_101_gf(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_101_gg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t2lecp

#endif /* LocalCorePotentialGeom101Func_hpp */
