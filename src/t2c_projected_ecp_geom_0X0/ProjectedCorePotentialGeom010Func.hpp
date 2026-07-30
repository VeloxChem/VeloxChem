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

#ifndef ProjectedCorePotentialGeom010Func_hpp
#define ProjectedCorePotentialGeom010Func_hpp

#include <array>
#include <iostream>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "ProjectedCorePotentialGeom010SSForS.hpp"
#include "ProjectedCorePotentialGeom010SSForP.hpp"
#include "ProjectedCorePotentialGeom010SSForD.hpp"
#include "ProjectedCorePotentialGeom010SSForF.hpp"
#include "ProjectedCorePotentialGeom010SSForG.hpp"
#include "ProjectedCorePotentialGeom010SPForS.hpp"
#include "ProjectedCorePotentialGeom010SPForP.hpp"
#include "ProjectedCorePotentialGeom010SPForD.hpp"
#include "ProjectedCorePotentialGeom010SPForF.hpp"
#include "ProjectedCorePotentialGeom010SPForG.hpp"
#include "ProjectedCorePotentialGeom010PSForS.hpp"
#include "ProjectedCorePotentialGeom010PSForP.hpp"
#include "ProjectedCorePotentialGeom010PSForD.hpp"
#include "ProjectedCorePotentialGeom010PSForF.hpp"
#include "ProjectedCorePotentialGeom010PSForG.hpp"
#include "ProjectedCorePotentialGeom010PPForS.hpp"
#include "ProjectedCorePotentialGeom010PPForP.hpp"
#include "ProjectedCorePotentialGeom010PPForD.hpp"
#include "ProjectedCorePotentialGeom010PPForF.hpp"
#include "ProjectedCorePotentialGeom010PPForG.hpp"

#include "ProjectedCorePotentialGeom010SDForS.hpp"
#include "ProjectedCorePotentialGeom010SDForP.hpp"
#include "ProjectedCorePotentialGeom010SDForD.hpp"
#include "ProjectedCorePotentialGeom010SDForF.hpp"
#include "ProjectedCorePotentialGeom010SDForG.hpp"
#include "ProjectedCorePotentialGeom010DSForS.hpp"
#include "ProjectedCorePotentialGeom010DSForP.hpp"
#include "ProjectedCorePotentialGeom010DSForD.hpp"
#include "ProjectedCorePotentialGeom010DSForF.hpp"
#include "ProjectedCorePotentialGeom010DSForG.hpp"

#include "ProjectedCorePotentialGeom010PDForS.hpp"
#include "ProjectedCorePotentialGeom010PDForP.hpp"
#include "ProjectedCorePotentialGeom010PDForD.hpp"
#include "ProjectedCorePotentialGeom010PDForF.hpp"
#include "ProjectedCorePotentialGeom010PDForG.hpp"
#include "ProjectedCorePotentialGeom010DPForS.hpp"
#include "ProjectedCorePotentialGeom010DPForP.hpp"
#include "ProjectedCorePotentialGeom010DPForD.hpp"
#include "ProjectedCorePotentialGeom010DPForF.hpp"
#include "ProjectedCorePotentialGeom010DPForG.hpp"

#include "ProjectedCorePotentialGeom010DDForS.hpp"
#include "ProjectedCorePotentialGeom010DDForP.hpp"
#include "ProjectedCorePotentialGeom010DDForD.hpp"
#include "ProjectedCorePotentialGeom010DDForF.hpp"
#include "ProjectedCorePotentialGeom010DDForG.hpp"

#include "ProjectedCorePotentialGeom010SFForS.hpp"
#include "ProjectedCorePotentialGeom010SFForP.hpp"
#include "ProjectedCorePotentialGeom010SFForD.hpp"
#include "ProjectedCorePotentialGeom010SFForF.hpp"
#include "ProjectedCorePotentialGeom010SFForG.hpp"
#include "ProjectedCorePotentialGeom010FSForS.hpp"
#include "ProjectedCorePotentialGeom010FSForP.hpp"
#include "ProjectedCorePotentialGeom010FSForD.hpp"
#include "ProjectedCorePotentialGeom010FSForF.hpp"
#include "ProjectedCorePotentialGeom010FSForG.hpp"

#include "ProjectedCorePotentialGeom010PFForS.hpp"
#include "ProjectedCorePotentialGeom010PFForP.hpp"
#include "ProjectedCorePotentialGeom010PFForD.hpp"
#include "ProjectedCorePotentialGeom010PFForF.hpp"
#include "ProjectedCorePotentialGeom010PFForG.hpp"
#include "ProjectedCorePotentialGeom010FPForS.hpp"
#include "ProjectedCorePotentialGeom010FPForP.hpp"
#include "ProjectedCorePotentialGeom010FPForD.hpp"
#include "ProjectedCorePotentialGeom010FPForF.hpp"
#include "ProjectedCorePotentialGeom010FPForG.hpp"

#include "ProjectedCorePotentialGeom010DFForS.hpp"
#include "ProjectedCorePotentialGeom010DFForP.hpp"
#include "ProjectedCorePotentialGeom010DFForD.hpp"
#include "ProjectedCorePotentialGeom010DFForF.hpp"
#include "ProjectedCorePotentialGeom010DFForG.hpp"
#include "ProjectedCorePotentialGeom010FDForS.hpp"
#include "ProjectedCorePotentialGeom010FDForP.hpp"
#include "ProjectedCorePotentialGeom010FDForD.hpp"
#include "ProjectedCorePotentialGeom010FDForF.hpp"
#include "ProjectedCorePotentialGeom010FDForG.hpp"

#include "ProjectedCorePotentialGeom010FFForS.hpp"
#include "ProjectedCorePotentialGeom010FFForP.hpp"
#include "ProjectedCorePotentialGeom010FFForD.hpp"
#include "ProjectedCorePotentialGeom010FFForF.hpp"
#include "ProjectedCorePotentialGeom010FFForG.hpp"

#include "ProjectedCorePotentialGeom010SGForS.hpp"
#include "ProjectedCorePotentialGeom010SGForP.hpp"
#include "ProjectedCorePotentialGeom010SGForD.hpp"
#include "ProjectedCorePotentialGeom010SGForF.hpp"
#include "ProjectedCorePotentialGeom010SGForG.hpp"
#include "ProjectedCorePotentialGeom010GSForS.hpp"
#include "ProjectedCorePotentialGeom010GSForP.hpp"
#include "ProjectedCorePotentialGeom010GSForD.hpp"
#include "ProjectedCorePotentialGeom010GSForF.hpp"
#include "ProjectedCorePotentialGeom010GSForG.hpp"

#include "ProjectedCorePotentialGeom010PGForS.hpp"
#include "ProjectedCorePotentialGeom010PGForP.hpp"
#include "ProjectedCorePotentialGeom010PGForD.hpp"
#include "ProjectedCorePotentialGeom010PGForF.hpp"
#include "ProjectedCorePotentialGeom010PGForG.hpp"
#include "ProjectedCorePotentialGeom010GPForS.hpp"
#include "ProjectedCorePotentialGeom010GPForP.hpp"
#include "ProjectedCorePotentialGeom010GPForD.hpp"
#include "ProjectedCorePotentialGeom010GPForF.hpp"
#include "ProjectedCorePotentialGeom010GPForG.hpp"

#include "ProjectedCorePotentialGeom010DGForS.hpp"
#include "ProjectedCorePotentialGeom010DGForP.hpp"
#include "ProjectedCorePotentialGeom010DGForD.hpp"
#include "ProjectedCorePotentialGeom010DGForF.hpp"
#include "ProjectedCorePotentialGeom010DGForG.hpp"
#include "ProjectedCorePotentialGeom010GDForS.hpp"
#include "ProjectedCorePotentialGeom010GDForP.hpp"
#include "ProjectedCorePotentialGeom010GDForD.hpp"
#include "ProjectedCorePotentialGeom010GDForF.hpp"
#include "ProjectedCorePotentialGeom010GDForG.hpp"

#include "ProjectedCorePotentialGeom010FGForS.hpp"
#include "ProjectedCorePotentialGeom010FGForP.hpp"
#include "ProjectedCorePotentialGeom010FGForD.hpp"
#include "ProjectedCorePotentialGeom010FGForF.hpp"
#include "ProjectedCorePotentialGeom010FGForG.hpp"
#include "ProjectedCorePotentialGeom010GFForS.hpp"
#include "ProjectedCorePotentialGeom010GFForP.hpp"
#include "ProjectedCorePotentialGeom010GFForD.hpp"
#include "ProjectedCorePotentialGeom010GFForF.hpp"
#include "ProjectedCorePotentialGeom010GFForG.hpp"

#include "ProjectedCorePotentialGeom010GGForS.hpp"
#include "ProjectedCorePotentialGeom010GGForP.hpp"
#include "ProjectedCorePotentialGeom010GGForD.hpp"
#include "ProjectedCorePotentialGeom010GGForF.hpp"
#include "ProjectedCorePotentialGeom010GGForG.hpp"

namespace t2pecp {

/// @brief Computes projected ECP integrals for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param ecp_potential The projected ECP potential.
/// @param ecp_momentum The momentum of projected ECP potential.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
compute_geom_010(T&                               distributor,
                 const CGtoBlock&                 bra_gto_block,
                 const CGtoBlock&                 ket_gto_block,
                 const CBaseCorePotential&        ecp_potential,
                 const int                        ecp_momentum,
                 const std::pair<size_t, size_t>& bra_indices,
                 const std::pair<size_t, size_t>& ket_indices,
                 const bool                       bra_eq_ket) -> void
{
    const auto bra_angmom = bra_gto_block.angular_momentum();

    const auto ket_angmom = ket_gto_block.angular_momentum();

    if ((bra_angmom == 0) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_ss_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_ss_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_ss_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_ss_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_ss_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_sp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_sp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_sp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_sp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_sp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_ps_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_ps_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_ps_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_ps_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_ps_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    
    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_pp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_pp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_pp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_pp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_pp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_sd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_sd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_sd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_sd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_sd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_ds_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_ds_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_ds_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_ds_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_ds_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_pd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_pd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_pd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_pd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_pd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_dp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_dp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_dp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_dp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_dp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_dd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_dd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_dd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_dd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_dd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_sf_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_sf_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_sf_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_sf_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_sf_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_fs_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_fs_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_fs_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_fs_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_fs_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_pf_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_pf_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_pf_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_pf_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_pf_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_fp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_fp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_fp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_fp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_fp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_df_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_df_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_df_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_df_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_df_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_fd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_fd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_fd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_fd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_fd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_ff_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_ff_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_ff_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_ff_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_ff_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_sg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_sg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_sg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_sg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_sg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_gs_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_gs_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_gs_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_gs_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_gs_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_pg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_pg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_pg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_pg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_pg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_gp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_gp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_gp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_gp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_gp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_dg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_dg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_dg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_dg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_dg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_gd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_gd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_gd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_gd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_gd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_fg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_fg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_fg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_fg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_fg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_gf_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_gf_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_gf_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_gf_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_gf_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_010_gg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_010_gg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_010_gg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_010_gg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_010_gg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
}

}  // namespace t2pecp

#endif /* ProjectedCorePotentialGeom010Func_hpp */
