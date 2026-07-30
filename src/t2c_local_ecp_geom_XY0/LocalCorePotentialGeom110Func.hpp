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

#ifndef LocalCorePotentialGeom110Func_hpp
#define LocalCorePotentialGeom110Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "LocalCorePotentialGeom110DD.hpp"
#include "LocalCorePotentialGeom110DF.hpp"
#include "LocalCorePotentialGeom110DG.hpp"
#include "LocalCorePotentialGeom110DP.hpp"
#include "LocalCorePotentialGeom110DS.hpp"
#include "LocalCorePotentialGeom110FD.hpp"
#include "LocalCorePotentialGeom110FF.hpp"
#include "LocalCorePotentialGeom110FG.hpp"
#include "LocalCorePotentialGeom110FP.hpp"
#include "LocalCorePotentialGeom110FS.hpp"
#include "LocalCorePotentialGeom110GD.hpp"
#include "LocalCorePotentialGeom110GF.hpp"
#include "LocalCorePotentialGeom110GG.hpp"
#include "LocalCorePotentialGeom110GP.hpp"
#include "LocalCorePotentialGeom110GS.hpp"
#include "LocalCorePotentialGeom110PD.hpp"
#include "LocalCorePotentialGeom110PF.hpp"
#include "LocalCorePotentialGeom110PG.hpp"
#include "LocalCorePotentialGeom110PP.hpp"
#include "LocalCorePotentialGeom110PS.hpp"
#include "LocalCorePotentialGeom110SD.hpp"
#include "LocalCorePotentialGeom110SF.hpp"
#include "LocalCorePotentialGeom110SG.hpp"
#include "LocalCorePotentialGeom110SP.hpp"
#include "LocalCorePotentialGeom110SS.hpp"

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
compute_geom_110(T&                               distributor,
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
        t2lecp::comp_local_core_potential_geom_110_ss(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_110_sp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        t2lecp::comp_local_core_potential_geom_110_ps(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_110_pp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_110_sd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        t2lecp::comp_local_core_potential_geom_110_ds(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_110_sf(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        t2lecp::comp_local_core_potential_geom_110_fs(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_110_pd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_110_dp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_110_sg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        t2lecp::comp_local_core_potential_geom_110_gs(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_110_pf(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_110_fp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_110_dd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_110_pg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_110_gp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_110_df(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_110_fd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_110_dg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_110_gd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_110_ff(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_110_fg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_110_gf(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_110_gg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t2lecp

#endif /* LocalCorePotentialGeom110Func_hpp */
