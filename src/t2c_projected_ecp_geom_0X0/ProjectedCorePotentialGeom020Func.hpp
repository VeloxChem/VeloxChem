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

#ifndef ProjectedCorePotentialGeom020Func_hpp
#define ProjectedCorePotentialGeom020Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "ProjectedCorePotentialGeom020SSForS.hpp"
#include "ProjectedCorePotentialGeom020SSForP.hpp"
#include "ProjectedCorePotentialGeom020SSForD.hpp"
#include "ProjectedCorePotentialGeom020SSForF.hpp"
#include "ProjectedCorePotentialGeom020SSForG.hpp"
#include "ProjectedCorePotentialGeom020SPForS.hpp"
#include "ProjectedCorePotentialGeom020SPForP.hpp"
#include "ProjectedCorePotentialGeom020SPForD.hpp"
#include "ProjectedCorePotentialGeom020SPForF.hpp"
#include "ProjectedCorePotentialGeom020SPForG.hpp"
#include "ProjectedCorePotentialGeom020PSForS.hpp"
#include "ProjectedCorePotentialGeom020PSForP.hpp"
#include "ProjectedCorePotentialGeom020PSForD.hpp"
#include "ProjectedCorePotentialGeom020PSForF.hpp"
#include "ProjectedCorePotentialGeom020PSForG.hpp"
#include "ProjectedCorePotentialGeom020PPForS.hpp"
#include "ProjectedCorePotentialGeom020PPForP.hpp"
#include "ProjectedCorePotentialGeom020PPForD.hpp"
#include "ProjectedCorePotentialGeom020PPForF.hpp"
#include "ProjectedCorePotentialGeom020PPForG.hpp"

#include "ProjectedCorePotentialGeom020SDForS.hpp"
#include "ProjectedCorePotentialGeom020SDForP.hpp"
#include "ProjectedCorePotentialGeom020SDForD.hpp"
#include "ProjectedCorePotentialGeom020SDForF.hpp"
#include "ProjectedCorePotentialGeom020SDForG.hpp"
#include "ProjectedCorePotentialGeom020DSForS.hpp"
#include "ProjectedCorePotentialGeom020DSForP.hpp"
#include "ProjectedCorePotentialGeom020DSForD.hpp"
#include "ProjectedCorePotentialGeom020DSForF.hpp"
#include "ProjectedCorePotentialGeom020DSForG.hpp"

#include "ProjectedCorePotentialGeom020PDForS.hpp"
#include "ProjectedCorePotentialGeom020PDForP.hpp"
#include "ProjectedCorePotentialGeom020PDForD.hpp"
#include "ProjectedCorePotentialGeom020PDForF.hpp"
#include "ProjectedCorePotentialGeom020PDForG.hpp"
#include "ProjectedCorePotentialGeom020DPForS.hpp"
#include "ProjectedCorePotentialGeom020DPForP.hpp"
#include "ProjectedCorePotentialGeom020DPForD.hpp"
#include "ProjectedCorePotentialGeom020DPForF.hpp"
#include "ProjectedCorePotentialGeom020DPForG.hpp"

#include "ProjectedCorePotentialGeom020DDForS.hpp"
#include "ProjectedCorePotentialGeom020DDForP.hpp"
#include "ProjectedCorePotentialGeom020DDForD.hpp"
#include "ProjectedCorePotentialGeom020DDForF.hpp"
#include "ProjectedCorePotentialGeom020DDForG.hpp"

#include "ProjectedCorePotentialGeom020SFForS.hpp"
#include "ProjectedCorePotentialGeom020SFForP.hpp"
#include "ProjectedCorePotentialGeom020SFForD.hpp"
#include "ProjectedCorePotentialGeom020SFForF.hpp"
#include "ProjectedCorePotentialGeom020SFForG.hpp"
#include "ProjectedCorePotentialGeom020FSForS.hpp"
#include "ProjectedCorePotentialGeom020FSForP.hpp"
#include "ProjectedCorePotentialGeom020FSForD.hpp"
#include "ProjectedCorePotentialGeom020FSForF.hpp"
#include "ProjectedCorePotentialGeom020FSForG.hpp"

#include "ProjectedCorePotentialGeom020PFForS.hpp"
#include "ProjectedCorePotentialGeom020PFForP.hpp"
#include "ProjectedCorePotentialGeom020PFForD.hpp"
#include "ProjectedCorePotentialGeom020PFForF.hpp"
#include "ProjectedCorePotentialGeom020PFForG.hpp"
#include "ProjectedCorePotentialGeom020FPForS.hpp"
#include "ProjectedCorePotentialGeom020FPForP.hpp"
#include "ProjectedCorePotentialGeom020FPForD.hpp"
#include "ProjectedCorePotentialGeom020FPForF.hpp"
#include "ProjectedCorePotentialGeom020FPForG.hpp"

#include "ProjectedCorePotentialGeom020DFForS.hpp"
#include "ProjectedCorePotentialGeom020DFForP.hpp"
#include "ProjectedCorePotentialGeom020DFForD.hpp"
#include "ProjectedCorePotentialGeom020DFForF.hpp"
#include "ProjectedCorePotentialGeom020DFForG.hpp"
#include "ProjectedCorePotentialGeom020FDForS.hpp"
#include "ProjectedCorePotentialGeom020FDForP.hpp"
#include "ProjectedCorePotentialGeom020FDForD.hpp"
#include "ProjectedCorePotentialGeom020FDForF.hpp"
#include "ProjectedCorePotentialGeom020FDForG.hpp"

#include "ProjectedCorePotentialGeom020FFForS.hpp"
#include "ProjectedCorePotentialGeom020FFForP.hpp"
#include "ProjectedCorePotentialGeom020FFForD.hpp"
#include "ProjectedCorePotentialGeom020FFForF.hpp"
#include "ProjectedCorePotentialGeom020FFForG.hpp"

#include "ProjectedCorePotentialGeom020SGForS.hpp"
#include "ProjectedCorePotentialGeom020SGForP.hpp"
#include "ProjectedCorePotentialGeom020SGForD.hpp"
#include "ProjectedCorePotentialGeom020SGForF.hpp"
#include "ProjectedCorePotentialGeom020SGForG.hpp"
#include "ProjectedCorePotentialGeom020GSForS.hpp"
#include "ProjectedCorePotentialGeom020GSForP.hpp"
#include "ProjectedCorePotentialGeom020GSForD.hpp"
#include "ProjectedCorePotentialGeom020GSForF.hpp"
#include "ProjectedCorePotentialGeom020GSForG.hpp"

#include "ProjectedCorePotentialGeom020PGForS.hpp"
#include "ProjectedCorePotentialGeom020PGForP.hpp"
#include "ProjectedCorePotentialGeom020PGForD.hpp"
#include "ProjectedCorePotentialGeom020PGForF.hpp"
#include "ProjectedCorePotentialGeom020PGForG.hpp"
#include "ProjectedCorePotentialGeom020GPForS.hpp"
#include "ProjectedCorePotentialGeom020GPForP.hpp"
#include "ProjectedCorePotentialGeom020GPForD.hpp"
#include "ProjectedCorePotentialGeom020GPForF.hpp"
#include "ProjectedCorePotentialGeom020GPForG.hpp"

#include "ProjectedCorePotentialGeom020DGForS.hpp"
#include "ProjectedCorePotentialGeom020DGForP.hpp"
#include "ProjectedCorePotentialGeom020DGForD.hpp"
#include "ProjectedCorePotentialGeom020DGForF.hpp"
#include "ProjectedCorePotentialGeom020DGForG.hpp"
#include "ProjectedCorePotentialGeom020GDForS.hpp"
#include "ProjectedCorePotentialGeom020GDForP.hpp"
#include "ProjectedCorePotentialGeom020GDForD.hpp"
#include "ProjectedCorePotentialGeom020GDForF.hpp"
#include "ProjectedCorePotentialGeom020GDForG.hpp"

#include "ProjectedCorePotentialGeom020FGForS.hpp"
#include "ProjectedCorePotentialGeom020FGForP.hpp"
#include "ProjectedCorePotentialGeom020FGForD.hpp"
#include "ProjectedCorePotentialGeom020FGForF.hpp"
#include "ProjectedCorePotentialGeom020FGForG.hpp"
#include "ProjectedCorePotentialGeom020GFForS.hpp"
#include "ProjectedCorePotentialGeom020GFForP.hpp"
#include "ProjectedCorePotentialGeom020GFForD.hpp"
#include "ProjectedCorePotentialGeom020GFForF.hpp"
#include "ProjectedCorePotentialGeom020GFForG.hpp"

#include "ProjectedCorePotentialGeom020GGForS.hpp"
#include "ProjectedCorePotentialGeom020GGForP.hpp"
#include "ProjectedCorePotentialGeom020GGForD.hpp"
#include "ProjectedCorePotentialGeom020GGForF.hpp"
#include "ProjectedCorePotentialGeom020GGForG.hpp"

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
compute_geom_020(T&                               distributor,
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
            t2pecp::comp_projected_core_potential_geom_020_ss_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_ss_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_ss_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_ss_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_ss_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_sp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_sp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_sp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_sp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_sp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_ps_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_ps_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_ps_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_ps_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_ps_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    
    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_pp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_pp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_pp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_pp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_pp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_sd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_sd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_sd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_sd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_sd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_ds_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_ds_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_ds_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_ds_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_ds_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_pd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_pd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_pd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_pd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_pd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_dp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_dp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_dp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_dp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_dp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_dd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_dd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_dd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_dd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_dd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_sf_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_sf_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_sf_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_sf_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_sf_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_fs_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_fs_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_fs_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_fs_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_fs_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_pf_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_pf_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_pf_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_pf_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_pf_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_fp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_fp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_fp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_fp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_fp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_df_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_df_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_df_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_df_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_df_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_fd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_fd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_fd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_fd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_fd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_ff_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_ff_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_ff_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_ff_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_ff_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_sg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_sg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_sg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_sg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_sg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_gs_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_gs_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_gs_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_gs_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_gs_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_pg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_pg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_pg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_pg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_pg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_gp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_gp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_gp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_gp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_gp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_dg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_dg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_dg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_dg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_dg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_gd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_gd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_gd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_gd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_gd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_fg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_fg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_fg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_fg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_fg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_gf_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_gf_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_gf_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_gf_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_gf_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_020_gg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_020_gg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_020_gg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_020_gg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_020_gg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
}

}  // namespace t2pecp

#endif /* ProjectedCorePotentialGeom020Func_hpp */
