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

#ifndef LocalCorePotentialGeom010Func_hpp
#define LocalCorePotentialGeom010Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "LocalCorePotentialGeom010DD.hpp"
#include "LocalCorePotentialGeom010DF.hpp"
#include "LocalCorePotentialGeom010DG.hpp"
#include "LocalCorePotentialGeom010DP.hpp"
#include "LocalCorePotentialGeom010DS.hpp"
#include "LocalCorePotentialGeom010FD.hpp"
#include "LocalCorePotentialGeom010FF.hpp"
#include "LocalCorePotentialGeom010FG.hpp"
#include "LocalCorePotentialGeom010FP.hpp"
#include "LocalCorePotentialGeom010FS.hpp"
#include "LocalCorePotentialGeom010GD.hpp"
#include "LocalCorePotentialGeom010GF.hpp"
#include "LocalCorePotentialGeom010GG.hpp"
#include "LocalCorePotentialGeom010GP.hpp"
#include "LocalCorePotentialGeom010GS.hpp"
#include "LocalCorePotentialGeom010PD.hpp"
#include "LocalCorePotentialGeom010PF.hpp"
#include "LocalCorePotentialGeom010PG.hpp"
#include "LocalCorePotentialGeom010PP.hpp"
#include "LocalCorePotentialGeom010PS.hpp"
#include "LocalCorePotentialGeom010SD.hpp"
#include "LocalCorePotentialGeom010SF.hpp"
#include "LocalCorePotentialGeom010SG.hpp"
#include "LocalCorePotentialGeom010SP.hpp"
#include "LocalCorePotentialGeom010SS.hpp"

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
compute_geom_010(T&                               distributor,
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
        t2lecp::comp_local_core_potential_geom_010_ss(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_010_sp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        t2lecp::comp_local_core_potential_geom_010_ps(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_010_pp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_010_sd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        t2lecp::comp_local_core_potential_geom_010_ds(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_010_sf(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        t2lecp::comp_local_core_potential_geom_010_fs(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_010_pd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_010_dp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_010_sg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        t2lecp::comp_local_core_potential_geom_010_gs(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_010_pf(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_010_fp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_010_dd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_010_pg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        t2lecp::comp_local_core_potential_geom_010_gp(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_010_df(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_010_fd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_010_dg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        t2lecp::comp_local_core_potential_geom_010_gd(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_010_ff(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_010_fg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        t2lecp::comp_local_core_potential_geom_010_gf(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        t2lecp::comp_local_core_potential_geom_010_gg(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t2lecp

#endif /* LocalCorePotentialGeom010Func_hpp */
