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




#include "SimdTwoCenterElectronRepulsionRsFunc.hpp"

#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"

#include "SimdElectronRepulsionRsRecSS.hpp"
#include "SimdElectronRepulsionRsRecSP.hpp"
#include "SimdElectronRepulsionRsRecSD.hpp"
#include "SimdElectronRepulsionRsRecSF.hpp"
#include "SimdElectronRepulsionRsRecSG.hpp"
#include "SimdElectronRepulsionRsRecSH.hpp"
#include "SimdElectronRepulsionRsRecSI.hpp"
#include "SimdElectronRepulsionRsRecSK.hpp"
#include "SimdElectronRepulsionRsRecSL.hpp"
#include "SimdElectronRepulsionRsRecPS.hpp"
#include "SimdElectronRepulsionRsRecPP.hpp"
#include "SimdElectronRepulsionRsRecPD.hpp"
#include "SimdElectronRepulsionRsRecPF.hpp"
#include "SimdElectronRepulsionRsRecPG.hpp"
#include "SimdElectronRepulsionRsRecPH.hpp"
#include "SimdElectronRepulsionRsRecPI.hpp"
#include "SimdElectronRepulsionRsRecPK.hpp"
#include "SimdElectronRepulsionRsRecPL.hpp"
#include "SimdElectronRepulsionRsRecDS.hpp"
#include "SimdElectronRepulsionRsRecDP.hpp"
#include "SimdElectronRepulsionRsRecDD.hpp"
#include "SimdElectronRepulsionRsRecDF.hpp"
#include "SimdElectronRepulsionRsRecDG.hpp"
#include "SimdElectronRepulsionRsRecDH.hpp"
#include "SimdElectronRepulsionRsRecDI.hpp"
#include "SimdElectronRepulsionRsRecDK.hpp"
#include "SimdElectronRepulsionRsRecDL.hpp"
#include "SimdElectronRepulsionRsRecFS.hpp"
#include "SimdElectronRepulsionRsRecFP.hpp"
#include "SimdElectronRepulsionRsRecFD.hpp"
#include "SimdElectronRepulsionRsRecFF.hpp"
#include "SimdElectronRepulsionRsRecFG.hpp"
#include "SimdElectronRepulsionRsRecFH.hpp"
#include "SimdElectronRepulsionRsRecFI.hpp"
#include "SimdElectronRepulsionRsRecFK.hpp"
#include "SimdElectronRepulsionRsRecFL.hpp"
#include "SimdElectronRepulsionRsRecGS.hpp"
#include "SimdElectronRepulsionRsRecGP.hpp"
#include "SimdElectronRepulsionRsRecGD.hpp"
#include "SimdElectronRepulsionRsRecGF.hpp"
#include "SimdElectronRepulsionRsRecGG.hpp"
#include "SimdElectronRepulsionRsRecGH.hpp"
#include "SimdElectronRepulsionRsRecGI.hpp"
#include "SimdElectronRepulsionRsRecGK.hpp"
#include "SimdElectronRepulsionRsRecGL.hpp"
#include "SimdElectronRepulsionRsRecHS.hpp"
#include "SimdElectronRepulsionRsRecHP.hpp"
#include "SimdElectronRepulsionRsRecHD.hpp"
#include "SimdElectronRepulsionRsRecHF.hpp"
#include "SimdElectronRepulsionRsRecHG.hpp"
#include "SimdElectronRepulsionRsRecHH.hpp"
#include "SimdElectronRepulsionRsRecHI.hpp"
#include "SimdElectronRepulsionRsRecHK.hpp"
#include "SimdElectronRepulsionRsRecHL.hpp"
#include "SimdElectronRepulsionRsRecIS.hpp"
#include "SimdElectronRepulsionRsRecIP.hpp"
#include "SimdElectronRepulsionRsRecID.hpp"
#include "SimdElectronRepulsionRsRecIF.hpp"
#include "SimdElectronRepulsionRsRecIG.hpp"
#include "SimdElectronRepulsionRsRecIH.hpp"
#include "SimdElectronRepulsionRsRecII.hpp"
#include "SimdElectronRepulsionRsRecIK.hpp"
#include "SimdElectronRepulsionRsRecIL.hpp"
#include "SimdElectronRepulsionRsRecKS.hpp"
#include "SimdElectronRepulsionRsRecKP.hpp"
#include "SimdElectronRepulsionRsRecKD.hpp"
#include "SimdElectronRepulsionRsRecKF.hpp"
#include "SimdElectronRepulsionRsRecKG.hpp"
#include "SimdElectronRepulsionRsRecKH.hpp"
#include "SimdElectronRepulsionRsRecKI.hpp"
#include "SimdElectronRepulsionRsRecKK.hpp"
#include "SimdElectronRepulsionRsRecKL.hpp"
#include "SimdElectronRepulsionRsRecLS.hpp"
#include "SimdElectronRepulsionRsRecLP.hpp"
#include "SimdElectronRepulsionRsRecLD.hpp"
#include "SimdElectronRepulsionRsRecLF.hpp"
#include "SimdElectronRepulsionRsRecLG.hpp"
#include "SimdElectronRepulsionRsRecLH.hpp"
#include "SimdElectronRepulsionRsRecLI.hpp"
#include "SimdElectronRepulsionRsRecLK.hpp"
#include "SimdElectronRepulsionRsRecLL.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
one_center_rs_electron_repulsion(const CBasisFunction &bra,
                                 const CBasisFunction &ket,
                                 const double          omega) -> std::pair<double, double>
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    if (lbra != lket) return {0.0, 0.0};

    // NOTE: the atoms meet, so the argument of the Boys function is zero and the
    // solid harmonic of the vector between them is one for the angular momentum
    // zero and vanishes above it. What is left is a closed formula in the
    // exponents alone, diagonal in the angular components and independent of
    // which component it is.

    constexpr auto fpi = mathconst::pi_value();

    const auto fcoul = 2.0 * fpi * fpi * std::sqrt(fpi);

    const auto fang = mathfunc::double_factorial(2 * lbra - 1) / static_cast<double>(2 * lbra + 1);

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    auto fsum = 0.0;

    auto fsum_rs = 0.0;

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

            const auto fval = a_norms[i] * b_norms[j] / fden;

            fsum += fval;

            // NOTE: the weight of the attenuation, which depends on the reduced
            // exponent of this pair of primitives and therefore belongs inside the
            // sum rather than outside it.

            const auto mu = a_exps[i] * b_exps[j] / fexp;

            const auto theta = omega / std::sqrt(omega * omega + mu);

            auto weight = theta;

            for (int l = 0; l < lbra; l++)
            {
                weight *= theta * theta;
            }

            fsum_rs += weight * fval;
        }
    }

    return {fcoul * fang * fsum, fcoul * fang * fsum_rs};
}

auto
compute_rs_electron_repulsion(double               *values,
                              const size_t          nvalues,
                              const CBasisFunction &bra,
                              const CBasisFunction &ket,
                              const CSimdMatrix    &coordinates,
                              CSimdMatrix          &buffer,
                              const double          omega) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    // NOTE: the combination is dispatched on the pair of angular momenta taken as a
    // single index, so the compiler forms one jump table rather than a chain of
    // comparisons which the combinations of high angular momentum reach last. The
    // unattenuated dispatcher does the same and for the same reason.

    switch (lbra * 9 + lket)
    {
        case   0:
            compute_rs_ss_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case   1:
            compute_rs_sp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case   2:
            compute_rs_sd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case   3:
            compute_rs_sf_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case   4:
            compute_rs_sg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case   5:
            compute_rs_sh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case   6:
            compute_rs_si_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case   7:
            compute_rs_sk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case   8:
            compute_rs_sl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case   9:
            compute_rs_ps_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  10:
            compute_rs_pp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  11:
            compute_rs_pd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  12:
            compute_rs_pf_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  13:
            compute_rs_pg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  14:
            compute_rs_ph_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  15:
            compute_rs_pi_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  16:
            compute_rs_pk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  17:
            compute_rs_pl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  18:
            compute_rs_ds_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  19:
            compute_rs_dp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  20:
            compute_rs_dd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  21:
            compute_rs_df_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  22:
            compute_rs_dg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  23:
            compute_rs_dh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  24:
            compute_rs_di_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  25:
            compute_rs_dk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  26:
            compute_rs_dl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  27:
            compute_rs_fs_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  28:
            compute_rs_fp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  29:
            compute_rs_fd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  30:
            compute_rs_ff_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  31:
            compute_rs_fg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  32:
            compute_rs_fh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  33:
            compute_rs_fi_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  34:
            compute_rs_fk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  35:
            compute_rs_fl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  36:
            compute_rs_gs_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  37:
            compute_rs_gp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  38:
            compute_rs_gd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  39:
            compute_rs_gf_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  40:
            compute_rs_gg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  41:
            compute_rs_gh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  42:
            compute_rs_gi_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  43:
            compute_rs_gk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  44:
            compute_rs_gl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  45:
            compute_rs_hs_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  46:
            compute_rs_hp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  47:
            compute_rs_hd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  48:
            compute_rs_hf_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  49:
            compute_rs_hg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  50:
            compute_rs_hh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  51:
            compute_rs_hi_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  52:
            compute_rs_hk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  53:
            compute_rs_hl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  54:
            compute_rs_is_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  55:
            compute_rs_ip_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  56:
            compute_rs_id_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  57:
            compute_rs_if_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  58:
            compute_rs_ig_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  59:
            compute_rs_ih_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  60:
            compute_rs_ii_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  61:
            compute_rs_ik_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  62:
            compute_rs_il_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  63:
            compute_rs_ks_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  64:
            compute_rs_kp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  65:
            compute_rs_kd_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  66:
            compute_rs_kf_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  67:
            compute_rs_kg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  68:
            compute_rs_kh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  69:
            compute_rs_ki_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  70:
            compute_rs_kk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  71:
            compute_rs_kl_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  72:
            compute_rs_ls_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  73:
            compute_rs_lp_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  74:
            compute_rs_ld_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  75:
            compute_rs_lf_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  76:
            compute_rs_lg_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  77:
            compute_rs_lh_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  78:
            compute_rs_li_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  79:
            compute_rs_lk_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        case  80:
            compute_rs_ll_electron_repulsion(values, nvalues, bra, ket, coordinates, buffer, omega);
            return;

        default:
            break;
    }

    // NOTE: the kernels are generated up to angular momentum eight, so a combination
    // above it stops rather than leaving the values of the matrix unwritten, which
    // is what a caller would otherwise read as integrals.

    errors::assertMsgCritical(
        false,
        std::string("SimdTwoCenterElectronRepulsionRsFunc.compute_rs_electron_repulsion: Integrals are not implemented"));
}

}  // namespace simdt2ceri
