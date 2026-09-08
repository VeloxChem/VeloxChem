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


#include "SimdElectronRepulsionRecFI.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdElectronRepulsionVrrRecDD.hpp"
#include "SimdElectronRepulsionVrrRecDF.hpp"
#include "SimdElectronRepulsionVrrRecDG.hpp"
#include "SimdElectronRepulsionVrrRecDH.hpp"
#include "SimdElectronRepulsionVrrRecDI.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFG.hpp"
#include "SimdElectronRepulsionVrrRecFH.hpp"
#include "SimdElectronRepulsionVrrRecFI.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPG.hpp"
#include "SimdElectronRepulsionVrrRecPH.hpp"
#include "SimdElectronRepulsionVrrRecPI.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSG.hpp"
#include "SimdElectronRepulsionVrrRecSH.hpp"
#include "SimdElectronRepulsionVrrRecSI.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdTransformFI.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_fi_electron_repulsion(double               *values,
                                   const size_t          nvalues,
                                   const CBasisFunction &bra,
                                   const CBasisFunction &ket,
                                   const CSimdMatrix    &coordinates) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_fi_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    auto buffer = CSimdMatrix(1725, nvalues);

    buffer.zero();

    const auto nmax = nvalues;

    for (size_t i = 0; i < nprim_a; i++)
    {
        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = nvalues;

            const auto p = a_exps[i] + b_exps[j];

            const auto mu = a_exps[i] * b_exps[j] / p;

            const auto pi = mathconst::pi_value();

            const auto fj = 2.0 * a_norms[i] * b_norms[j] * pi * pi * std::sqrt(pi)
                            / (a_exps[i] * b_exps[j] * std::sqrt(p));

            const auto alpha = a_exps[i];

            const auto beta = b_exps[j];

            const auto fb = a_exps[i] / p;

            const auto fa = -b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            simdfunc::compute_pb(buffer, coordinates, 3, ncols, fb);

            simdfunc::compute_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7, 8, 9}, ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 16, 0, 7, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 19, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 14, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 15, ncols);

            compute_prim_ds_electron_repulsion_1(buffer, 43, 0, 7, 8, 22, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 46, 0, 8, 9, 25, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 49, 0, 9, 10, 28, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 52, 0, 10, 11, 31, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 55, 0, 11, 12, 34, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 58, 0, 12, 13, 37, ncols, alpha, beta, p);

            compute_prim_ds_electron_repulsion_1(buffer, 61, 0, 13, 14, 40, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 64, 0, 16, 19, 43, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 67, 0, 19, 22, 46, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 70, 0, 22, 25, 49, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 73, 0, 25, 28, 52, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 76, 0, 28, 31, 55, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 79, 0, 31, 34, 58, ncols, alpha, beta, p);

            compute_prim_fs_electron_repulsion_1(buffer, 82, 0, 34, 37, 61, ncols, alpha, beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 85, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 88, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 91, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 94, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 97, 3, 14, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 100, 3, 15, ncols);

            compute_prim_pp_electron_repulsion_2(buffer, 103, 3, 9, 25, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 106, 3, 10, 28, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 109, 3, 11, 31, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 112, 3, 12, 34, ncols, p);

            compute_prim_pp_electron_repulsion_2(buffer, 115, 3, 13, 37, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 118, 3, 22, 46, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 127, 3, 25, 49, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 136, 3, 28, 52, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 145, 3, 31, 55, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 154, 3, 34, 58, ncols, p);

            compute_prim_dp_electron_repulsion_2(buffer, 163, 3, 37, 61, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 172, 3, 46, 70, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 181, 3, 49, 73, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 190, 3, 52, 76, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 199, 3, 55, 79, ncols, p);

            compute_prim_fp_electron_repulsion_2(buffer, 208, 3, 58, 82, ncols, p);

            compute_prim_sd_electron_repulsion_1(buffer, 217, 3, 9, 10, 88, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 220, 3, 10, 11, 91, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 223, 3, 11, 12, 94, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 226, 3, 12, 13, 97, ncols, alpha, beta, p);

            compute_prim_sd_electron_repulsion_1(buffer, 229, 3, 13, 14, 100, ncols, alpha, beta, p);

            compute_prim_pd_electron_repulsion_4(buffer, 232, 0, 85, 217, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 235, 0, 88, 220, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 238, 0, 91, 223, ncols, p);

            compute_prim_pd_electron_repulsion_4(buffer, 241, 0, 94, 226, ncols, p);

            compute_prim_dd_electron_repulsion_2(buffer, 244, 3, 103, 43, 46, 127, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 253, 3, 106, 46, 49, 136, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 262, 3, 109, 49, 52, 145, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 271, 3, 112, 52, 55, 154, ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_2(buffer, 280, 3, 115, 55, 58, 163, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 289, 3, 118, 64, 67, 172, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 298, 3, 127, 67, 70, 181, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 307, 3, 136, 70, 73, 190, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 316, 3, 145, 73, 76, 199, ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_2(buffer, 325, 3, 154, 76, 79, 208, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_3(buffer, 334, 3, 85, 88, 220, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_3(buffer, 343, 3, 88, 91, 223, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 352, 3, 91, 94, 226, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_1(buffer, 358, 3, 94, 97, 229, ncols, alpha, beta, p);

            compute_prim_pf_electron_repulsion_9(buffer, 364, 0, 217, 334, ncols, p);

            compute_prim_pf_electron_repulsion_13(buffer, 367, 0, 3, 220, 343, 238, ncols, p);

            compute_prim_pf_electron_repulsion_8(buffer, 371, 0, 223, 352, ncols, p);

            compute_prim_df_electron_repulsion_2(buffer, 374, 3, 232, 118, 127, 253, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 392, 3, 235, 127, 136, 262, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 410, 3, 238, 136, 145, 271, ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_2(buffer, 428, 3, 241, 145, 154, 280, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 446, 3, 253, 172, 181, 307, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 464, 3, 262, 181, 190, 316, ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_2(buffer, 482, 3, 271, 190, 199, 325, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_3(buffer, 500, 3, 217, 220, 343, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_4(buffer, 512, 3, 220, 223, 352, ncols, alpha, beta, p);

            compute_prim_sg_electron_repulsion_1(buffer, 524, 3, 223, 226, 358, ncols, alpha, beta, p);

            compute_prim_pg_electron_repulsion_11(buffer, 533, 0, 3, 334, 500, 367, ncols, p);

            compute_prim_pg_electron_repulsion_12(buffer, 545, 0, 3, 343, 512, 371, ncols, p);

            compute_prim_pg_electron_repulsion_13(buffer, 557, 0, 352, 524, ncols, p);

            compute_prim_dg_electron_repulsion_2(buffer, 560, 3, 364, 244, 253, 392, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_7(buffer, 587, 0, 3, 367, 545, 253, 262, 410, ncols, alpha, beta, p);

            compute_prim_dg_electron_repulsion_8(buffer, 624, 0, 3, 371, 557, 262, 271, 428, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_2(buffer, 654, 3, 374, 289, 298, 446, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_2(buffer, 681, 3, 392, 298, 307, 464, ncols, alpha, beta, p);

            compute_prim_fg_electron_repulsion_3(buffer, 708, 0, 3, 410, 624, 307, 316, 482, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_7(buffer, 750, 3, 334, 343, 512, ncols, alpha, beta, p);

            compute_prim_sh_electron_repulsion_6(buffer, 763, 3, 343, 352, 524, ncols, alpha, beta, p);

            compute_prim_ph_electron_repulsion_6(buffer, 776, 0, 3, 500, 750, 545, ncols, p);

            compute_prim_ph_electron_repulsion_7(buffer, 804, 0, 512, 763, ncols, p);

            compute_prim_dh_electron_repulsion_3(buffer, 813, 0, 3, 533, 776, 374, 392, 587, ncols, alpha, beta, p);

            compute_prim_dh_electron_repulsion_4(buffer, 892, 0, 3, 545, 804, 392, 410, 624, ncols, alpha, beta, p);

            compute_prim_fh_electron_repulsion_1(buffer, 949, 0, 3, 587, 892, 446, 464, 708, ncols, alpha, beta, p);

            compute_prim_si_electron_repulsion_1(buffer, 1075, 3, 500, 512, 763, ncols, alpha, beta, p);

            compute_prim_pi_electron_repulsion_2(buffer, 1088, 0, 750, 1075, ncols, p);

            compute_prim_di_electron_repulsion_1(buffer, 1101, 0, 3, 776, 1088, 560, 587, 892, ncols, alpha, beta, p);

            compute_prim_fi_electron_repulsion_0(buffer, 1165, 0, 3, 813, 1101, 654, 681, 949, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1445, 1165, 280, ncols);
        }
    }

    simdtrf::transform_fi(values, nvalues, buffer, 1445, nmax);
}

}  // namespace simdt2ceri
