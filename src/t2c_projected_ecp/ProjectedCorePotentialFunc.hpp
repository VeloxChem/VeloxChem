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

#ifndef ProjectedCorePotentialFunc_hpp
#define ProjectedCorePotentialFunc_hpp

#include <array>

#include "GtoBlock.hpp"
#include "BaseCorePotential.hpp"
#include "ProjectedCorePotentialSSForS.hpp"
#include "ProjectedCorePotentialSSForP.hpp"
#include "ProjectedCorePotentialSSForD.hpp"
#include "ProjectedCorePotentialSSForF.hpp"
#include "ProjectedCorePotentialSSForG.hpp"
#include "ProjectedCorePotentialSPForS.hpp"
#include "ProjectedCorePotentialSPForP.hpp"
#include "ProjectedCorePotentialSPForD.hpp"

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
compute(T&                               distributor,
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
            t2pecp::comp_projected_core_potential_ss_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_ss_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_ss_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 3)
        {
            t2pecp::comp_projected_core_potential_ss_for_f(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 4)
        {
            t2pecp::comp_projected_core_potential_ss_for_g(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        if (ecp_momentum == 0)
        {
            t2pecp::comp_projected_core_potential_sp_for_s(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 1)
        {
            t2pecp::comp_projected_core_potential_sp_for_p(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        if (ecp_momentum == 2)
        {
            t2pecp::comp_projected_core_potential_sp_for_d(distributor, bra_gto_block, ket_gto_block, ecp_potential, bra_indices, ket_indices, bra_eq_ket);

            return;
        }
        
        std::cout << " *** ECP projectors of angular momentum " << ecp_momentum << " are not supported. " << std::endl;
    }
}

}  // namespace t2lecp

#endif /* ProjectedCorePotentialFunc_hpp */
