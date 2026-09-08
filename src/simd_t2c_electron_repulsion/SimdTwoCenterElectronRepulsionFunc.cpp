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



#include "SimdTwoCenterElectronRepulsionFunc.hpp"

#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "SimdElectronRepulsionRecDD.hpp"
#include "SimdElectronRepulsionRecDF.hpp"
#include "SimdElectronRepulsionRecDG.hpp"
#include "SimdElectronRepulsionRecDH.hpp"
#include "SimdElectronRepulsionRecDI.hpp"
#include "SimdElectronRepulsionRecDP.hpp"
#include "SimdElectronRepulsionRecDS.hpp"
#include "SimdElectronRepulsionRecFD.hpp"
#include "SimdElectronRepulsionRecFF.hpp"
#include "SimdElectronRepulsionRecFG.hpp"
#include "SimdElectronRepulsionRecFH.hpp"
#include "SimdElectronRepulsionRecFI.hpp"
#include "SimdElectronRepulsionRecFP.hpp"
#include "SimdElectronRepulsionRecFS.hpp"
#include "SimdElectronRepulsionRecGD.hpp"
#include "SimdElectronRepulsionRecGF.hpp"
#include "SimdElectronRepulsionRecGG.hpp"
#include "SimdElectronRepulsionRecGH.hpp"
#include "SimdElectronRepulsionRecGI.hpp"
#include "SimdElectronRepulsionRecGP.hpp"
#include "SimdElectronRepulsionRecGS.hpp"
#include "SimdElectronRepulsionRecHD.hpp"
#include "SimdElectronRepulsionRecHF.hpp"
#include "SimdElectronRepulsionRecHG.hpp"
#include "SimdElectronRepulsionRecHH.hpp"
#include "SimdElectronRepulsionRecHI.hpp"
#include "SimdElectronRepulsionRecHP.hpp"
#include "SimdElectronRepulsionRecHS.hpp"
#include "SimdElectronRepulsionRecID.hpp"
#include "SimdElectronRepulsionRecIF.hpp"
#include "SimdElectronRepulsionRecIG.hpp"
#include "SimdElectronRepulsionRecIH.hpp"
#include "SimdElectronRepulsionRecII.hpp"
#include "SimdElectronRepulsionRecIP.hpp"
#include "SimdElectronRepulsionRecIS.hpp"
#include "SimdElectronRepulsionRecPD.hpp"
#include "SimdElectronRepulsionRecPF.hpp"
#include "SimdElectronRepulsionRecPG.hpp"
#include "SimdElectronRepulsionRecPH.hpp"
#include "SimdElectronRepulsionRecPI.hpp"
#include "SimdElectronRepulsionRecPP.hpp"
#include "SimdElectronRepulsionRecPS.hpp"
#include "SimdElectronRepulsionRecSD.hpp"
#include "SimdElectronRepulsionRecSF.hpp"
#include "SimdElectronRepulsionRecSG.hpp"
#include "SimdElectronRepulsionRecSH.hpp"
#include "SimdElectronRepulsionRecSI.hpp"
#include "SimdElectronRepulsionRecSP.hpp"
#include "SimdElectronRepulsionRecSS.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
one_center_electron_repulsion(const CBasisFunction &bra, const CBasisFunction &ket) -> double
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    if (lbra != lket) return 0.0;

    // NOTE: the atoms meet, so the argument of the Boys function is zero and the
    // solid harmonic of the vector between them is one for the angular momentum
    // zero and vanishes above it. What is left is a closed formula in the
    // exponents alone, diagonal in the angular components and independent of
    // which component it is.

    constexpr auto fpi = mathconst::pi_value();

    const auto fcoul = 2.0 * fpi * fpi * std::sqrt(fpi);

    // NOTE: the double factorial of twice the angular momentum less one over
    // twice it plus one is the factor the angular momentum contributes, and is
    // one for the angular momentum zero.

    const auto fang = mathfunc::double_factorial(2 * lbra - 1) / static_cast<double>(2 * lbra + 1);

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    auto fsum = 0.0;

    for (size_t i = 0; i < a_exps.size(); i++)
    {
        for (size_t j = 0; j < b_exps.size(); j++)
        {
            const auto fexp = a_exps[i] + b_exps[j];

            auto fden = a_exps[i] * b_exps[j] * std::sqrt(fexp);

            // NOTE: twice the sum of the exponents is raised to the angular
            // momentum by repeated multiplication, as the angular momentum is
            // small.

            for (int l = 0; l < lbra; l++)
            {
                fden *= 2.0 * fexp;
            }

            fsum += a_norms[i] * b_norms[j] / fden;
        }
    }

    return fcoul * fang * fsum;
}

auto
compute_electron_repulsion(double               *values,
                           const size_t          nvalues,
                           const CBasisFunction &bra,
                           const CBasisFunction &ket,
                           const CSimdMatrix    &coordinates) -> void
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
            compute_ss_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case  1:
            compute_sp_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case  2:
            compute_sd_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case  3:
            compute_sf_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case  4:
            compute_sg_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case  5:
            compute_sh_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case  6:
            compute_si_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case  7:
            compute_ps_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case  8:
            compute_pp_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case  9:
            compute_pd_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 10:
            compute_pf_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 11:
            compute_pg_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 12:
            compute_ph_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 13:
            compute_pi_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 14:
            compute_ds_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 15:
            compute_dp_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 16:
            compute_dd_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 17:
            compute_df_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 18:
            compute_dg_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 19:
            compute_dh_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 20:
            compute_di_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 21:
            compute_fs_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 22:
            compute_fp_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 23:
            compute_fd_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 24:
            compute_ff_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 25:
            compute_fg_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 26:
            compute_fh_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 27:
            compute_fi_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 28:
            compute_gs_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 29:
            compute_gp_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 30:
            compute_gd_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 31:
            compute_gf_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 32:
            compute_gg_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 33:
            compute_gh_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 34:
            compute_gi_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 35:
            compute_hs_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 36:
            compute_hp_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 37:
            compute_hd_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 38:
            compute_hf_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 39:
            compute_hg_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 40:
            compute_hh_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 41:
            compute_hi_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 42:
            compute_is_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 43:
            compute_ip_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 44:
            compute_id_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 45:
            compute_if_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 46:
            compute_ig_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 47:
            compute_ih_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        case 48:
            compute_ii_electron_repulsion(values, nvalues, bra, ket, coordinates);
            return;

        default:
            break;
    }

    // NOTE: the kernels are generated up to angular momentum six, so a combination
    // above it stops rather than leaving the values of the matrix unwritten, which
    // is what a caller would otherwise read as integrals.

    errors::assertMsgCritical(
        false, std::string("SimdTwoCenterElectronRepulsionFunc.compute_electron_repulsion: Integrals are not implemented"));
}

}  // namespace simdt2ceri
