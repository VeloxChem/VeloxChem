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


#include "SimdNuclearPotentialFunc.hpp"

#include <string>

#include "ErrorHandler.hpp"

#include "SimdNuclearPotentialRecSS.hpp"
#include "SimdNuclearPotentialRecSP.hpp"
#include "SimdNuclearPotentialRecSD.hpp"
#include "SimdNuclearPotentialRecSF.hpp"
#include "SimdNuclearPotentialRecSG.hpp"
#include "SimdNuclearPotentialRecSH.hpp"
#include "SimdNuclearPotentialRecSI.hpp"
#include "SimdNuclearPotentialRecPS.hpp"
#include "SimdNuclearPotentialRecPP.hpp"
#include "SimdNuclearPotentialRecPD.hpp"
#include "SimdNuclearPotentialRecPF.hpp"
#include "SimdNuclearPotentialRecPG.hpp"
#include "SimdNuclearPotentialRecPH.hpp"
#include "SimdNuclearPotentialRecPI.hpp"
#include "SimdNuclearPotentialRecDS.hpp"
#include "SimdNuclearPotentialRecDP.hpp"
#include "SimdNuclearPotentialRecDD.hpp"
#include "SimdNuclearPotentialRecDF.hpp"
#include "SimdNuclearPotentialRecDG.hpp"
#include "SimdNuclearPotentialRecDH.hpp"
#include "SimdNuclearPotentialRecDI.hpp"
#include "SimdNuclearPotentialRecFS.hpp"
#include "SimdNuclearPotentialRecFP.hpp"
#include "SimdNuclearPotentialRecFD.hpp"
#include "SimdNuclearPotentialRecFF.hpp"
#include "SimdNuclearPotentialRecFG.hpp"
#include "SimdNuclearPotentialRecFH.hpp"
#include "SimdNuclearPotentialRecFI.hpp"
#include "SimdNuclearPotentialRecGS.hpp"
#include "SimdNuclearPotentialRecGP.hpp"
#include "SimdNuclearPotentialRecGD.hpp"
#include "SimdNuclearPotentialRecGF.hpp"
#include "SimdNuclearPotentialRecGG.hpp"
#include "SimdNuclearPotentialRecGH.hpp"
#include "SimdNuclearPotentialRecGI.hpp"
#include "SimdNuclearPotentialRecHS.hpp"
#include "SimdNuclearPotentialRecHP.hpp"
#include "SimdNuclearPotentialRecHD.hpp"
#include "SimdNuclearPotentialRecHF.hpp"
#include "SimdNuclearPotentialRecHG.hpp"
#include "SimdNuclearPotentialRecHH.hpp"
#include "SimdNuclearPotentialRecHI.hpp"
#include "SimdNuclearPotentialRecIS.hpp"
#include "SimdNuclearPotentialRecIP.hpp"
#include "SimdNuclearPotentialRecID.hpp"
#include "SimdNuclearPotentialRecIF.hpp"
#include "SimdNuclearPotentialRecIG.hpp"
#include "SimdNuclearPotentialRecIH.hpp"
#include "SimdNuclearPotentialRecII.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_nuclear_potential(double                    *values,
                          const size_t               nvalues,
                          const CBasisFunction      &bra,
                          const CBasisFunction      &ket,
                          const CSimdMatrix         &coordinates,
                          const std::vector<double> &charges,
                          const std::vector<double> &points,
                          CSimdMatrix               &buffer,
                          const double               threshold) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    // NOTE: the two orders of a combination are separate kernels, as the recurrence
    // builds the angular momentum on ket side and transfers it to bra side, so the
    // work differs with the order even where the integrals do not.

    // NOTE: the combination is dispatched on the pair of angular momenta taken as a
    // single index, so the compiler forms one jump table rather than a chain of
    // comparisons which the combinations of high angular momentum reach last.

    switch (lbra * 7 + lket)
    {
        case   0:
            compute_ss_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case   1:
            compute_sp_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case   2:
            compute_sd_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case   3:
            compute_sf_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case   4:
            compute_sg_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case   5:
            compute_sh_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case   6:
            compute_si_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case   7:
            compute_ps_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case   8:
            compute_pp_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case   9:
            compute_pd_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  10:
            compute_pf_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  11:
            compute_pg_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  12:
            compute_ph_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  13:
            compute_pi_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  14:
            compute_ds_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  15:
            compute_dp_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  16:
            compute_dd_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  17:
            compute_df_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  18:
            compute_dg_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  19:
            compute_dh_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  20:
            compute_di_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  21:
            compute_fs_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  22:
            compute_fp_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  23:
            compute_fd_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  24:
            compute_ff_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  25:
            compute_fg_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  26:
            compute_fh_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  27:
            compute_fi_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  28:
            compute_gs_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  29:
            compute_gp_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  30:
            compute_gd_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  31:
            compute_gf_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  32:
            compute_gg_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  33:
            compute_gh_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  34:
            compute_gi_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  35:
            compute_hs_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  36:
            compute_hp_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  37:
            compute_hd_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  38:
            compute_hf_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  39:
            compute_hg_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  40:
            compute_hh_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  41:
            compute_hi_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  42:
            compute_is_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  43:
            compute_ip_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  44:
            compute_id_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  45:
            compute_if_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  46:
            compute_ig_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  47:
            compute_ih_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;

        case  48:
            compute_ii_nuclear_potential(values, nvalues, bra, ket, coordinates, charges, points, buffer, threshold);
            return;
        default:
            break;
    }

    // NOTE: falling through to nothing would leave the values as the caller found
    // them and say nothing, and a caller cannot tell integrals which were not
    // computed from integrals which are zero.

    errors::assertMsgCritical(false,
                              std::string("compute_nuclear_potential: No kernel for the combination of angular momenta ") +
                                  std::to_string(lbra) + std::string(" and ") + std::to_string(lket) +
                                  std::string("; the kernels reach angular momentum six"));
}

}  // namespace simdnpot
