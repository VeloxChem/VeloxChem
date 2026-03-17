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

#ifndef ProjectedCorePotentialGeom110Func_hpp
#define ProjectedCorePotentialGeom110Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "ProjectedCorePotentialGeom110SSForS.hpp"
#include "ProjectedCorePotentialGeom110SSForP.hpp"
#include "ProjectedCorePotentialGeom110SSForD.hpp"
#include "ProjectedCorePotentialGeom110SSForF.hpp"
#include "ProjectedCorePotentialGeom110SSForG.hpp"
#include "ProjectedCorePotentialGeom110SPForS.hpp"
#include "ProjectedCorePotentialGeom110SPForP.hpp"
#include "ProjectedCorePotentialGeom110SPForD.hpp"
#include "ProjectedCorePotentialGeom110SPForF.hpp"
#include "ProjectedCorePotentialGeom110SPForG.hpp"
#include "ProjectedCorePotentialGeom110PSForS.hpp"
#include "ProjectedCorePotentialGeom110PSForP.hpp"
#include "ProjectedCorePotentialGeom110PSForD.hpp"
#include "ProjectedCorePotentialGeom110PSForF.hpp"
#include "ProjectedCorePotentialGeom110PSForG.hpp"
#include "ProjectedCorePotentialGeom110PPForS.hpp"
#include "ProjectedCorePotentialGeom110PPForP.hpp"
#include "ProjectedCorePotentialGeom110PPForD.hpp"
#include "ProjectedCorePotentialGeom110PPForF.hpp"
#include "ProjectedCorePotentialGeom110PPForG.hpp"

#include "ProjectedCorePotentialGeom110SDForS.hpp"
#include "ProjectedCorePotentialGeom110SDForP.hpp"
#include "ProjectedCorePotentialGeom110SDForD.hpp"
#include "ProjectedCorePotentialGeom110SDForF.hpp"
#include "ProjectedCorePotentialGeom110SDForG.hpp"
#include "ProjectedCorePotentialGeom110DSForS.hpp"
#include "ProjectedCorePotentialGeom110DSForP.hpp"
#include "ProjectedCorePotentialGeom110DSForD.hpp"
#include "ProjectedCorePotentialGeom110DSForF.hpp"
#include "ProjectedCorePotentialGeom110DSForG.hpp"

#include "ProjectedCorePotentialGeom110PDForS.hpp"
#include "ProjectedCorePotentialGeom110PDForP.hpp"
#include "ProjectedCorePotentialGeom110PDForD.hpp"
#include "ProjectedCorePotentialGeom110PDForF.hpp"
#include "ProjectedCorePotentialGeom110PDForG.hpp"
#include "ProjectedCorePotentialGeom110DPForS.hpp"
#include "ProjectedCorePotentialGeom110DPForP.hpp"
#include "ProjectedCorePotentialGeom110DPForD.hpp"
#include "ProjectedCorePotentialGeom110DPForF.hpp"
#include "ProjectedCorePotentialGeom110DPForG.hpp"

#include "ProjectedCorePotentialGeom110DDForS.hpp"
#include "ProjectedCorePotentialGeom110DDForP.hpp"
#include "ProjectedCorePotentialGeom110DDForD.hpp"
#include "ProjectedCorePotentialGeom110DDForF.hpp"
#include "ProjectedCorePotentialGeom110DDForG.hpp"

#include "ProjectedCorePotentialGeom110SFForS.hpp"
#include "ProjectedCorePotentialGeom110SFForP.hpp"
#include "ProjectedCorePotentialGeom110SFForD.hpp"
#include "ProjectedCorePotentialGeom110SFForF.hpp"
#include "ProjectedCorePotentialGeom110SFForG.hpp"
#include "ProjectedCorePotentialGeom110FSForS.hpp"
#include "ProjectedCorePotentialGeom110FSForP.hpp"
#include "ProjectedCorePotentialGeom110FSForD.hpp"
#include "ProjectedCorePotentialGeom110FSForF.hpp"
#include "ProjectedCorePotentialGeom110FSForG.hpp"

#include "ProjectedCorePotentialGeom110PFForS.hpp"
#include "ProjectedCorePotentialGeom110PFForP.hpp"
#include "ProjectedCorePotentialGeom110PFForD.hpp"
#include "ProjectedCorePotentialGeom110PFForF.hpp"
#include "ProjectedCorePotentialGeom110PFForG.hpp"
#include "ProjectedCorePotentialGeom110FPForS.hpp"
#include "ProjectedCorePotentialGeom110FPForP.hpp"
#include "ProjectedCorePotentialGeom110FPForD.hpp"
#include "ProjectedCorePotentialGeom110FPForF.hpp"
#include "ProjectedCorePotentialGeom110FPForG.hpp"

#include "ProjectedCorePotentialGeom110DFForS.hpp"
#include "ProjectedCorePotentialGeom110DFForP.hpp"
#include "ProjectedCorePotentialGeom110DFForD.hpp"
#include "ProjectedCorePotentialGeom110DFForF.hpp"
#include "ProjectedCorePotentialGeom110DFForG.hpp"
#include "ProjectedCorePotentialGeom110FDForS.hpp"
#include "ProjectedCorePotentialGeom110FDForP.hpp"
#include "ProjectedCorePotentialGeom110FDForD.hpp"
#include "ProjectedCorePotentialGeom110FDForF.hpp"
#include "ProjectedCorePotentialGeom110FDForG.hpp"

#include "ProjectedCorePotentialGeom110FFForS.hpp"
#include "ProjectedCorePotentialGeom110FFForP.hpp"
#include "ProjectedCorePotentialGeom110FFForD.hpp"
#include "ProjectedCorePotentialGeom110FFForF.hpp"
#include "ProjectedCorePotentialGeom110FFForG.hpp"

#include "ProjectedCorePotentialGeom110SGForS.hpp"
#include "ProjectedCorePotentialGeom110SGForP.hpp"
#include "ProjectedCorePotentialGeom110SGForD.hpp"
#include "ProjectedCorePotentialGeom110SGForF.hpp"
#include "ProjectedCorePotentialGeom110SGForG.hpp"
#include "ProjectedCorePotentialGeom110GSForS.hpp"
#include "ProjectedCorePotentialGeom110GSForP.hpp"
#include "ProjectedCorePotentialGeom110GSForD.hpp"
#include "ProjectedCorePotentialGeom110GSForF.hpp"
#include "ProjectedCorePotentialGeom110GSForG.hpp"

#include "ProjectedCorePotentialGeom110PGForS.hpp"
#include "ProjectedCorePotentialGeom110PGForP.hpp"
#include "ProjectedCorePotentialGeom110PGForD.hpp"
#include "ProjectedCorePotentialGeom110PGForF.hpp"
#include "ProjectedCorePotentialGeom110PGForG.hpp"
#include "ProjectedCorePotentialGeom110GPForS.hpp"
#include "ProjectedCorePotentialGeom110GPForP.hpp"
#include "ProjectedCorePotentialGeom110GPForD.hpp"
#include "ProjectedCorePotentialGeom110GPForF.hpp"
#include "ProjectedCorePotentialGeom110GPForG.hpp"

#include "ProjectedCorePotentialGeom110DGForS.hpp"
#include "ProjectedCorePotentialGeom110DGForP.hpp"
#include "ProjectedCorePotentialGeom110DGForD.hpp"
#include "ProjectedCorePotentialGeom110DGForF.hpp"
#include "ProjectedCorePotentialGeom110DGForG.hpp"
#include "ProjectedCorePotentialGeom110GDForS.hpp"
#include "ProjectedCorePotentialGeom110GDForP.hpp"
#include "ProjectedCorePotentialGeom110GDForD.hpp"
#include "ProjectedCorePotentialGeom110GDForF.hpp"
#include "ProjectedCorePotentialGeom110GDForG.hpp"

#include "ProjectedCorePotentialGeom110FGForS.hpp"
#include "ProjectedCorePotentialGeom110FGForP.hpp"
#include "ProjectedCorePotentialGeom110FGForD.hpp"
#include "ProjectedCorePotentialGeom110FGForF.hpp"
#include "ProjectedCorePotentialGeom110FGForG.hpp"
#include "ProjectedCorePotentialGeom110GFForS.hpp"
#include "ProjectedCorePotentialGeom110GFForP.hpp"
#include "ProjectedCorePotentialGeom110GFForD.hpp"
#include "ProjectedCorePotentialGeom110GFForF.hpp"
#include "ProjectedCorePotentialGeom110GFForG.hpp"

#include "ProjectedCorePotentialGeom110GGForS.hpp"
#include "ProjectedCorePotentialGeom110GGForP.hpp"
#include "ProjectedCorePotentialGeom110GGForD.hpp"
#include "ProjectedCorePotentialGeom110GGForF.hpp"
#include "ProjectedCorePotentialGeom110GGForG.hpp"

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
compute_geom_110(T&                               distributor,
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
            t2pecp::comp_projected_core_potential_geom_110_ss_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_ss_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_ss_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_ss_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_ss_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_sp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_sp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_sp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_sp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_sp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_ps_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_ps_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_ps_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_ps_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_ps_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    
    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_pp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_pp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_pp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_pp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_pp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_sd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_sd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_sd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_sd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_sd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_ds_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_ds_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_ds_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_ds_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_ds_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_pd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_pd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_pd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_pd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_pd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_dp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_dp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_dp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_dp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_dp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_dd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_dd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_dd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_dd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_dd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_sf_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_sf_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_sf_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_sf_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_sf_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_fs_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_fs_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_fs_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_fs_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_fs_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_pf_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_pf_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_pf_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_pf_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_pf_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_fp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_fp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_fp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_fp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_fp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_df_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_df_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_df_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_df_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_df_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_fd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_fd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_fd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_fd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_fd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_ff_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_ff_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_ff_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_ff_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_ff_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_sg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_sg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_sg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_sg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_sg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_gs_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_gs_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_gs_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_gs_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_gs_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_pg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_pg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_pg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_pg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_pg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_gp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_gp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_gp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_gp_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_gp_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_dg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_dg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_dg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_dg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_dg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_gd_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_gd_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_gd_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_gd_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_gd_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_fg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_fg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_fg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_fg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_fg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_gf_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_gf_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_gf_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_gf_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_gf_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_geom_110_gg_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_geom_110_gg_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_geom_110_gg_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_geom_110_gg_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_geom_110_gg_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
}

}  // namespace t2pecp

#endif /* ProjectedCorePotentialGeom110Func_hpp */
