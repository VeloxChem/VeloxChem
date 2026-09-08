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



#include "SimdKineticEnergyFunc.hpp"

#include <string>

#include "ErrorHandler.hpp"
#include "SimdKineticEnergyRecDD.hpp"
#include "SimdKineticEnergyRecDF.hpp"
#include "SimdKineticEnergyRecDG.hpp"
#include "SimdKineticEnergyRecDH.hpp"
#include "SimdKineticEnergyRecDI.hpp"
#include "SimdKineticEnergyRecDP.hpp"
#include "SimdKineticEnergyRecDS.hpp"
#include "SimdKineticEnergyRecFD.hpp"
#include "SimdKineticEnergyRecFF.hpp"
#include "SimdKineticEnergyRecFG.hpp"
#include "SimdKineticEnergyRecFH.hpp"
#include "SimdKineticEnergyRecFI.hpp"
#include "SimdKineticEnergyRecFP.hpp"
#include "SimdKineticEnergyRecFS.hpp"
#include "SimdKineticEnergyRecGD.hpp"
#include "SimdKineticEnergyRecGF.hpp"
#include "SimdKineticEnergyRecGG.hpp"
#include "SimdKineticEnergyRecGH.hpp"
#include "SimdKineticEnergyRecGI.hpp"
#include "SimdKineticEnergyRecGP.hpp"
#include "SimdKineticEnergyRecGS.hpp"
#include "SimdKineticEnergyRecHD.hpp"
#include "SimdKineticEnergyRecHF.hpp"
#include "SimdKineticEnergyRecHG.hpp"
#include "SimdKineticEnergyRecHH.hpp"
#include "SimdKineticEnergyRecHI.hpp"
#include "SimdKineticEnergyRecHP.hpp"
#include "SimdKineticEnergyRecHS.hpp"
#include "SimdKineticEnergyRecID.hpp"
#include "SimdKineticEnergyRecIF.hpp"
#include "SimdKineticEnergyRecIG.hpp"
#include "SimdKineticEnergyRecIH.hpp"
#include "SimdKineticEnergyRecII.hpp"
#include "SimdKineticEnergyRecIP.hpp"
#include "SimdKineticEnergyRecIS.hpp"
#include "SimdKineticEnergyRecPD.hpp"
#include "SimdKineticEnergyRecPF.hpp"
#include "SimdKineticEnergyRecPG.hpp"
#include "SimdKineticEnergyRecPH.hpp"
#include "SimdKineticEnergyRecPI.hpp"
#include "SimdKineticEnergyRecPP.hpp"
#include "SimdKineticEnergyRecPS.hpp"
#include "SimdKineticEnergyRecSD.hpp"
#include "SimdKineticEnergyRecSF.hpp"
#include "SimdKineticEnergyRecSG.hpp"
#include "SimdKineticEnergyRecSH.hpp"
#include "SimdKineticEnergyRecSI.hpp"
#include "SimdKineticEnergyRecSP.hpp"
#include "SimdKineticEnergyRecSS.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_kinetic_energy(double               *values,
                       const size_t          nvalues,
                       const CBasisFunction &bra,
                       const CBasisFunction &ket,
                       const CSimdMatrix    &coordinates,
                       const double          threshold) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    // NOTE: the two orders of a combination are separate kernels, as the recurrence
    // builds the angular momentum on one side and reaches the other through it, so
    // the work differs with the order even where the integrals do not.

    // NOTE: the combination is dispatched on the pair of angular momenta taken as a
    // single index, so the compiler forms one jump table rather than a chain of
    // comparisons which the combinations of high angular momentum reach last.

    switch (lbra * 7 + lket)
    {
        case  0:
            compute_ss_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  1:
            compute_sp_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  2:
            compute_sd_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  3:
            compute_sf_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  4:
            compute_sg_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  5:
            compute_sh_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  6:
            compute_si_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  7:
            compute_ps_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  8:
            compute_pp_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  9:
            compute_pd_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 10:
            compute_pf_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 11:
            compute_pg_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 12:
            compute_ph_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 13:
            compute_pi_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 14:
            compute_ds_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 15:
            compute_dp_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 16:
            compute_dd_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 17:
            compute_df_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 18:
            compute_dg_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 19:
            compute_dh_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 20:
            compute_di_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 21:
            compute_fs_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 22:
            compute_fp_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 23:
            compute_fd_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 24:
            compute_ff_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 25:
            compute_fg_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 26:
            compute_fh_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 27:
            compute_fi_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 28:
            compute_gs_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 29:
            compute_gp_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 30:
            compute_gd_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 31:
            compute_gf_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 32:
            compute_gg_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 33:
            compute_gh_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 34:
            compute_gi_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 35:
            compute_hs_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 36:
            compute_hp_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 37:
            compute_hd_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 38:
            compute_hf_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 39:
            compute_hg_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 40:
            compute_hh_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 41:
            compute_hi_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 42:
            compute_is_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 43:
            compute_ip_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 44:
            compute_id_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 45:
            compute_if_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 46:
            compute_ig_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 47:
            compute_ih_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 48:
            compute_ii_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);
            return;

        default:
            break;
    }

    // NOTE: the kernels are generated up to angular momentum six, so a combination
    // above it stops rather than leaving the values of the sparsity pattern
    // unwritten, which is what a caller would otherwise read as integrals.

    errors::assertMsgCritical(false, std::string("SimdKineticEnergyFunc.compute_kinetic_energy: Kinetic energy integrals are not implemented"));
}

}  // namespace simdkin
