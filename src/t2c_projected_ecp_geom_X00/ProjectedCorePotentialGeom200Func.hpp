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

#ifndef ProjectedCorePotentialGeom200Func_hpp
#define ProjectedCorePotentialGeom200Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "ProjectedCorePotentialGeom200SSForS.hpp"
#include "ProjectedCorePotentialGeom200SSForP.hpp"
#include "ProjectedCorePotentialGeom200SSForD.hpp"
#include "ProjectedCorePotentialGeom200SSForF.hpp"
#include "ProjectedCorePotentialGeom200SSForG.hpp"
#include "ProjectedCorePotentialGeom200SPForS.hpp"
#include "ProjectedCorePotentialGeom200SPForP.hpp"
#include "ProjectedCorePotentialGeom200SPForD.hpp"
#include "ProjectedCorePotentialGeom200SPForF.hpp"
#include "ProjectedCorePotentialGeom200SPForG.hpp"
#include "ProjectedCorePotentialGeom200PSForS.hpp"
#include "ProjectedCorePotentialGeom200PSForP.hpp"
#include "ProjectedCorePotentialGeom200PSForD.hpp"
#include "ProjectedCorePotentialGeom200PSForF.hpp"
#include "ProjectedCorePotentialGeom200PSForG.hpp"
#include "ProjectedCorePotentialGeom200PPForS.hpp"
#include "ProjectedCorePotentialGeom200PPForP.hpp"
#include "ProjectedCorePotentialGeom200PPForD.hpp"
#include "ProjectedCorePotentialGeom200PPForF.hpp"
#include "ProjectedCorePotentialGeom200PPForG.hpp"

#include "ProjectedCorePotentialGeom200SDForS.hpp"
#include "ProjectedCorePotentialGeom200SDForP.hpp"
#include "ProjectedCorePotentialGeom200SDForD.hpp"
#include "ProjectedCorePotentialGeom200SDForF.hpp"
#include "ProjectedCorePotentialGeom200SDForG.hpp"
#include "ProjectedCorePotentialGeom200DSForS.hpp"
#include "ProjectedCorePotentialGeom200DSForP.hpp"
#include "ProjectedCorePotentialGeom200DSForD.hpp"
#include "ProjectedCorePotentialGeom200DSForF.hpp"
#include "ProjectedCorePotentialGeom200DSForG.hpp"

#include "ProjectedCorePotentialGeom200PDForS.hpp"
#include "ProjectedCorePotentialGeom200PDForP.hpp"
#include "ProjectedCorePotentialGeom200PDForD.hpp"
#include "ProjectedCorePotentialGeom200PDForF.hpp"
#include "ProjectedCorePotentialGeom200PDForG.hpp"
#include "ProjectedCorePotentialGeom200DPForS.hpp"
#include "ProjectedCorePotentialGeom200DPForP.hpp"
#include "ProjectedCorePotentialGeom200DPForD.hpp"
#include "ProjectedCorePotentialGeom200DPForF.hpp"
#include "ProjectedCorePotentialGeom200DPForG.hpp"

#include "ProjectedCorePotentialGeom200DDForS.hpp"
#include "ProjectedCorePotentialGeom200DDForP.hpp"
#include "ProjectedCorePotentialGeom200DDForD.hpp"
#include "ProjectedCorePotentialGeom200DDForF.hpp"
#include "ProjectedCorePotentialGeom200DDForG.hpp"

#include "ProjectedCorePotentialGeom200SFForS.hpp"
#include "ProjectedCorePotentialGeom200SFForP.hpp"
#include "ProjectedCorePotentialGeom200SFForD.hpp"
#include "ProjectedCorePotentialGeom200SFForF.hpp"
#include "ProjectedCorePotentialGeom200SFForG.hpp"
#include "ProjectedCorePotentialGeom200FSForS.hpp"
#include "ProjectedCorePotentialGeom200FSForP.hpp"
#include "ProjectedCorePotentialGeom200FSForD.hpp"
#include "ProjectedCorePotentialGeom200FSForF.hpp"
#include "ProjectedCorePotentialGeom200FSForG.hpp"

#include "ProjectedCorePotentialGeom200PFForS.hpp"
#include "ProjectedCorePotentialGeom200PFForP.hpp"
#include "ProjectedCorePotentialGeom200PFForD.hpp"
#include "ProjectedCorePotentialGeom200PFForF.hpp"
#include "ProjectedCorePotentialGeom200PFForG.hpp"
#include "ProjectedCorePotentialGeom200FPForS.hpp"
#include "ProjectedCorePotentialGeom200FPForP.hpp"
#include "ProjectedCorePotentialGeom200FPForD.hpp"
#include "ProjectedCorePotentialGeom200FPForF.hpp"
#include "ProjectedCorePotentialGeom200FPForG.hpp"

#include "ProjectedCorePotentialGeom200DFForS.hpp"
#include "ProjectedCorePotentialGeom200DFForP.hpp"
#include "ProjectedCorePotentialGeom200DFForD.hpp"
#include "ProjectedCorePotentialGeom200DFForF.hpp"
#include "ProjectedCorePotentialGeom200DFForG.hpp"
#include "ProjectedCorePotentialGeom200FDForS.hpp"
#include "ProjectedCorePotentialGeom200FDForP.hpp"
#include "ProjectedCorePotentialGeom200FDForD.hpp"
#include "ProjectedCorePotentialGeom200FDForF.hpp"
#include "ProjectedCorePotentialGeom200FDForG.hpp"

#include "ProjectedCorePotentialGeom200FFForS.hpp"
#include "ProjectedCorePotentialGeom200FFForP.hpp"
#include "ProjectedCorePotentialGeom200FFForD.hpp"
#include "ProjectedCorePotentialGeom200FFForF.hpp"
#include "ProjectedCorePotentialGeom200FFForG.hpp"

#include "ProjectedCorePotentialGeom200SGForS.hpp"
#include "ProjectedCorePotentialGeom200SGForP.hpp"
#include "ProjectedCorePotentialGeom200SGForD.hpp"
#include "ProjectedCorePotentialGeom200SGForF.hpp"
#include "ProjectedCorePotentialGeom200SGForG.hpp"
#include "ProjectedCorePotentialGeom200GSForS.hpp"
#include "ProjectedCorePotentialGeom200GSForP.hpp"
#include "ProjectedCorePotentialGeom200GSForD.hpp"
#include "ProjectedCorePotentialGeom200GSForF.hpp"
#include "ProjectedCorePotentialGeom200GSForG.hpp"

#include "ProjectedCorePotentialGeom200PGForS.hpp"
#include "ProjectedCorePotentialGeom200PGForP.hpp"
#include "ProjectedCorePotentialGeom200PGForD.hpp"
#include "ProjectedCorePotentialGeom200PGForF.hpp"
#include "ProjectedCorePotentialGeom200PGForG.hpp"
#include "ProjectedCorePotentialGeom200GPForS.hpp"
#include "ProjectedCorePotentialGeom200GPForP.hpp"
#include "ProjectedCorePotentialGeom200GPForD.hpp"
#include "ProjectedCorePotentialGeom200GPForF.hpp"
#include "ProjectedCorePotentialGeom200GPForG.hpp"

#include "ProjectedCorePotentialGeom200DGForS.hpp"
#include "ProjectedCorePotentialGeom200DGForP.hpp"
#include "ProjectedCorePotentialGeom200DGForD.hpp"
#include "ProjectedCorePotentialGeom200DGForF.hpp"
#include "ProjectedCorePotentialGeom200DGForG.hpp"
#include "ProjectedCorePotentialGeom200GDForS.hpp"
#include "ProjectedCorePotentialGeom200GDForP.hpp"
#include "ProjectedCorePotentialGeom200GDForD.hpp"
#include "ProjectedCorePotentialGeom200GDForF.hpp"
#include "ProjectedCorePotentialGeom200GDForG.hpp"

#include "ProjectedCorePotentialGeom200FGForS.hpp"
#include "ProjectedCorePotentialGeom200FGForP.hpp"
#include "ProjectedCorePotentialGeom200FGForD.hpp"
#include "ProjectedCorePotentialGeom200FGForF.hpp"
#include "ProjectedCorePotentialGeom200FGForG.hpp"
#include "ProjectedCorePotentialGeom200GFForS.hpp"
#include "ProjectedCorePotentialGeom200GFForP.hpp"
#include "ProjectedCorePotentialGeom200GFForD.hpp"
#include "ProjectedCorePotentialGeom200GFForF.hpp"
#include "ProjectedCorePotentialGeom200GFForG.hpp"

#include "ProjectedCorePotentialGeom200GGForS.hpp"
#include "ProjectedCorePotentialGeom200GGForP.hpp"
#include "ProjectedCorePotentialGeom200GGForD.hpp"
#include "ProjectedCorePotentialGeom200GGForF.hpp"
#include "ProjectedCorePotentialGeom200GGForG.hpp"

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
compute_geom_200(T&                               distributor,
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
            t2pecp::comp_projected_core_potential_geom_200_ss_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_ss_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_ss_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_ss_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_ss_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_sp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_sp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_sp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_sp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_sp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_ps_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_ps_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_ps_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_ps_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_ps_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    
    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_pp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_pp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_pp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_pp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_pp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_sd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_sd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_sd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_sd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_sd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_ds_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_ds_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_ds_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_ds_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_ds_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_pd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_pd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_pd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_pd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_pd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_dp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_dp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_dp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_dp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_dp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_dd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_dd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_dd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_dd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_dd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_sf_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_sf_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_sf_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_sf_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_sf_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_fs_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_fs_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_fs_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_fs_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_fs_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_pf_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_pf_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_pf_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_pf_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_pf_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_fp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_fp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_fp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_fp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_fp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_df_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_df_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_df_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_df_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_df_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_fd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_fd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_fd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_fd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_fd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_ff_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_ff_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_ff_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_ff_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_ff_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_sg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_sg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_sg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_sg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_sg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_gs_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_gs_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_gs_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_gs_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_gs_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_pg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_pg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_pg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_pg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_pg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_gp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_gp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_gp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_gp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_gp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_dg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_dg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_dg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_dg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_dg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_gd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_gd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_gd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_gd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_gd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_fg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_fg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_fg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_fg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_fg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_gf_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_gf_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_gf_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_gf_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_gf_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_200_gg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_200_gg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_200_gg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_200_gg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_200_gg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
}

}  // namespace t2pecp

#endif /* ProjectedCorePotentialGeom200Func_hpp */
