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

#ifndef DiagonalElectronRepulsionFunc_hpp
#define DiagonalElectronRepulsionFunc_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionDiagRecDDDD.hpp"
#include "ElectronRepulsionDiagRecDFDF.hpp"
#include "ElectronRepulsionDiagRecDGDG.hpp"
#include "ElectronRepulsionDiagRecFFFF.hpp"
#include "ElectronRepulsionDiagRecFGFG.hpp"
#include "ElectronRepulsionDiagRecGGGG.hpp"
#include "ElectronRepulsionDiagRecPDPD.hpp"
#include "ElectronRepulsionDiagRecPFPF.hpp"
#include "ElectronRepulsionDiagRecPGPG.hpp"
#include "ElectronRepulsionDiagRecPPPP.hpp"
#include "ElectronRepulsionDiagRecSDSD.hpp"
#include "ElectronRepulsionDiagRecSFSF.hpp"
#include "ElectronRepulsionDiagRecSGSG.hpp"
#include "ElectronRepulsionDiagRecSPSP.hpp"
#include "ElectronRepulsionDiagRecSSSS.hpp"
#include "GtoPairBlock.hpp"

namespace erifunc {  // erifunc namespace

/// @brief Computes vector of largest of Cartesian integral components of requested four center integral.
/// @param distributor  The pointer to distributor of Fock matrix or Fock matrices.
/// @param gto_pair_block  The GTOs pair block to compute Fock matrix or Fock matrices.
/// @param range  The range [first, last) for bra and ket computations.
template <class T>
inline auto
diag_compute(T& distributor, const CGtoPairBlock& gto_pair_block, const std::pair<size_t, size_t>& range) -> void
{
    const auto angmoms = gto_pair_block.angular_momentums();

    // angular order - 0

    if (angmoms == std::pair<int, int>({0, 0}))
    {
        erirec::comp_diag_electron_repulsion_ssss(distributor, gto_pair_block, range);

        return;
    }

    // angular order - 1

    if (angmoms == std::pair<int, int>({0, 1}))
    {
        erirec::comp_diag_electron_repulsion_spsp(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({1, 0}))
    {
        erirec::comp_diag_electron_repulsion_psps(distributor, gto_pair_block, range);

        return;
    }

    // angular order - 2

    if (angmoms == std::pair<int, int>({0, 2}))
    {
        erirec::comp_diag_electron_repulsion_sdsd(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({2, 0}))
    {
        erirec::comp_diag_electron_repulsion_dsds(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({1, 1}))
    {
        erirec::comp_diag_electron_repulsion_pppp(distributor, gto_pair_block, range);

        return;
    }

    // angular order - 3

    if (angmoms == std::pair<int, int>({0, 3}))
    {
        erirec::comp_diag_electron_repulsion_sfsf(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({3, 0}))
    {
        erirec::comp_diag_electron_repulsion_fsfs(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({1, 2}))
    {
        erirec::comp_diag_electron_repulsion_pdpd(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({2, 1}))
    {
        erirec::comp_diag_electron_repulsion_dpdp(distributor, gto_pair_block, range);

        return;
    }

    // angular order - 4

    if (angmoms == std::pair<int, int>({0, 4}))
    {
        erirec::comp_diag_electron_repulsion_sgsg(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({4, 0}))
    {
        erirec::comp_diag_electron_repulsion_gsgs(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({1, 3}))
    {
        erirec::comp_diag_electron_repulsion_pfpf(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({3, 1}))
    {
        erirec::comp_diag_electron_repulsion_fpfp(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({2, 2}))
    {
        erirec::comp_diag_electron_repulsion_dddd(distributor, gto_pair_block, range);

        return;
    }

    // angular order - 5

    if (angmoms == std::pair<int, int>({1, 4}))
    {
        erirec::comp_diag_electron_repulsion_pgpg(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({4, 1}))
    {
        erirec::comp_diag_electron_repulsion_gpgp(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({2, 3}))
    {
        erirec::comp_diag_electron_repulsion_dfdf(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({3, 2}))
    {
        erirec::comp_diag_electron_repulsion_fdfd(distributor, gto_pair_block, range);

        return;
    }

    // angular order - 6

    if (angmoms == std::pair<int, int>({2, 4}))
    {
        erirec::comp_diag_electron_repulsion_dgdg(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({4, 2}))
    {
        erirec::comp_diag_electron_repulsion_gdgd(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({3, 3}))
    {
        erirec::comp_diag_electron_repulsion_ffff(distributor, gto_pair_block, range);

        return;
    }

    // angular order - 7

    if (angmoms == std::pair<int, int>({3, 4}))
    {
        erirec::comp_diag_electron_repulsion_fgfg(distributor, gto_pair_block, range);

        return;
    }

    if (angmoms == std::pair<int, int>({3, 4}))
    {
        erirec::comp_diag_electron_repulsion_gfgf(distributor, gto_pair_block, range);

        return;
    }

    // angular order - 8

    if (angmoms == std::pair<int, int>({4, 4}))
    {
        erirec::comp_diag_electron_repulsion_gggg(distributor, gto_pair_block, range);

        return;
    }
}

}  // namespace erifunc

#endif /* DiagonalElectronRepulsionFunc_hpp */
