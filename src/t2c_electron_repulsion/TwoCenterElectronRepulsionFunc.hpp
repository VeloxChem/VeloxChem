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

#ifndef TwoCenterElectronRepulsionFunc_hpp
#define TwoCenterElectronRepulsionFunc_hpp

#include <array>

#include "GtoBlock.hpp"
#include "TwoCenterElectronRepulsionRecII.hpp"
#include "TwoCenterElectronRepulsionRecIH.hpp"
#include "TwoCenterElectronRepulsionRecID.hpp"
#include "TwoCenterElectronRepulsionRecIF.hpp"
#include "TwoCenterElectronRepulsionRecIG.hpp"
#include "TwoCenterElectronRepulsionRecIP.hpp"
#include "TwoCenterElectronRepulsionRecIS.hpp"
#include "TwoCenterElectronRepulsionRecHI.hpp"
#include "TwoCenterElectronRepulsionRecHH.hpp"
#include "TwoCenterElectronRepulsionRecHD.hpp"
#include "TwoCenterElectronRepulsionRecHF.hpp"
#include "TwoCenterElectronRepulsionRecHG.hpp"
#include "TwoCenterElectronRepulsionRecHP.hpp"
#include "TwoCenterElectronRepulsionRecHS.hpp"
#include "TwoCenterElectronRepulsionRecDI.hpp"
#include "TwoCenterElectronRepulsionRecDH.hpp"
#include "TwoCenterElectronRepulsionRecDD.hpp"
#include "TwoCenterElectronRepulsionRecDF.hpp"
#include "TwoCenterElectronRepulsionRecDG.hpp"
#include "TwoCenterElectronRepulsionRecDP.hpp"
#include "TwoCenterElectronRepulsionRecDS.hpp"
#include "TwoCenterElectronRepulsionRecFI.hpp"
#include "TwoCenterElectronRepulsionRecFH.hpp"
#include "TwoCenterElectronRepulsionRecFD.hpp"
#include "TwoCenterElectronRepulsionRecFF.hpp"
#include "TwoCenterElectronRepulsionRecFG.hpp"
#include "TwoCenterElectronRepulsionRecFP.hpp"
#include "TwoCenterElectronRepulsionRecFS.hpp"
#include "TwoCenterElectronRepulsionRecGI.hpp"
#include "TwoCenterElectronRepulsionRecGH.hpp"
#include "TwoCenterElectronRepulsionRecGD.hpp"
#include "TwoCenterElectronRepulsionRecGF.hpp"
#include "TwoCenterElectronRepulsionRecGG.hpp"
#include "TwoCenterElectronRepulsionRecGP.hpp"
#include "TwoCenterElectronRepulsionRecGS.hpp"
#include "TwoCenterElectronRepulsionRecPI.hpp"
#include "TwoCenterElectronRepulsionRecPH.hpp"
#include "TwoCenterElectronRepulsionRecPD.hpp"
#include "TwoCenterElectronRepulsionRecPF.hpp"
#include "TwoCenterElectronRepulsionRecPG.hpp"
#include "TwoCenterElectronRepulsionRecPP.hpp"
#include "TwoCenterElectronRepulsionRecPS.hpp"
#include "TwoCenterElectronRepulsionRecSI.hpp"
#include "TwoCenterElectronRepulsionRecSH.hpp"
#include "TwoCenterElectronRepulsionRecSD.hpp"
#include "TwoCenterElectronRepulsionRecSF.hpp"
#include "TwoCenterElectronRepulsionRecSG.hpp"
#include "TwoCenterElectronRepulsionRecSP.hpp"
#include "TwoCenterElectronRepulsionRecSS.hpp"

namespace t2cerifunc {

/// @brief Computes electron repulsion integrals for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
compute(T&                               distributor,
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
        t2ceri::comp_electron_repulsion_ss(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_sp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_ps(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_pp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_sd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_ds(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_sf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_fs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_pd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_dp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_sg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_gs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_pf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_fp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_dd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_sh(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 5) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_hs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_pg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_gp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_df(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_fd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_si(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 6) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_is(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_ph(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 5) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_hp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_dg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_gd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_ff(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_pi(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 6) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_ip(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_dh(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 5) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_hd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_fg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_gf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_di(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 6) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_id(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_fh(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 5) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_hf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_gg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_fi(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 6) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_if(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_gh(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 5) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_hg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_gi(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 6) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_ig(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 5) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_hh(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 5) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_hi(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 6) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_ih(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 6) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_ii(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t2cerifunc

#endif /* TwoCenterElectronRepulsionFunc_hpp */
