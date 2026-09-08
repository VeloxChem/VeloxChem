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


#include "SimdKineticEnergyRecHI.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdKineticEnergyVrrRecDD.hpp"
#include "SimdKineticEnergyVrrRecDF.hpp"
#include "SimdKineticEnergyVrrRecDG.hpp"
#include "SimdKineticEnergyVrrRecDH.hpp"
#include "SimdKineticEnergyVrrRecDI.hpp"
#include "SimdKineticEnergyVrrRecDP.hpp"
#include "SimdKineticEnergyVrrRecDS.hpp"
#include "SimdKineticEnergyVrrRecFD.hpp"
#include "SimdKineticEnergyVrrRecFF.hpp"
#include "SimdKineticEnergyVrrRecFG.hpp"
#include "SimdKineticEnergyVrrRecFH.hpp"
#include "SimdKineticEnergyVrrRecFI.hpp"
#include "SimdKineticEnergyVrrRecFP.hpp"
#include "SimdKineticEnergyVrrRecFS.hpp"
#include "SimdKineticEnergyVrrRecGD.hpp"
#include "SimdKineticEnergyVrrRecGF.hpp"
#include "SimdKineticEnergyVrrRecGG.hpp"
#include "SimdKineticEnergyVrrRecGH.hpp"
#include "SimdKineticEnergyVrrRecGI.hpp"
#include "SimdKineticEnergyVrrRecGP.hpp"
#include "SimdKineticEnergyVrrRecGS.hpp"
#include "SimdKineticEnergyVrrRecHD.hpp"
#include "SimdKineticEnergyVrrRecHF.hpp"
#include "SimdKineticEnergyVrrRecHG.hpp"
#include "SimdKineticEnergyVrrRecHH.hpp"
#include "SimdKineticEnergyVrrRecHI.hpp"
#include "SimdKineticEnergyVrrRecHP.hpp"
#include "SimdKineticEnergyVrrRecHS.hpp"
#include "SimdKineticEnergyVrrRecPD.hpp"
#include "SimdKineticEnergyVrrRecPF.hpp"
#include "SimdKineticEnergyVrrRecPG.hpp"
#include "SimdKineticEnergyVrrRecPH.hpp"
#include "SimdKineticEnergyVrrRecPI.hpp"
#include "SimdKineticEnergyVrrRecPP.hpp"
#include "SimdKineticEnergyVrrRecPS.hpp"
#include "SimdKineticEnergyVrrRecSD.hpp"
#include "SimdKineticEnergyVrrRecSF.hpp"
#include "SimdKineticEnergyVrrRecSG.hpp"
#include "SimdKineticEnergyVrrRecSH.hpp"
#include "SimdKineticEnergyVrrRecSI.hpp"
#include "SimdKineticEnergyVrrRecSP.hpp"
#include "SimdKineticEnergyVrrRecSS.hpp"
#include "SimdOverlapVrrRecDD.hpp"
#include "SimdOverlapVrrRecDF.hpp"
#include "SimdOverlapVrrRecDG.hpp"
#include "SimdOverlapVrrRecDH.hpp"
#include "SimdOverlapVrrRecDI.hpp"
#include "SimdOverlapVrrRecDP.hpp"
#include "SimdOverlapVrrRecDS.hpp"
#include "SimdOverlapVrrRecFD.hpp"
#include "SimdOverlapVrrRecFF.hpp"
#include "SimdOverlapVrrRecFG.hpp"
#include "SimdOverlapVrrRecFH.hpp"
#include "SimdOverlapVrrRecFI.hpp"
#include "SimdOverlapVrrRecFP.hpp"
#include "SimdOverlapVrrRecFS.hpp"
#include "SimdOverlapVrrRecGD.hpp"
#include "SimdOverlapVrrRecGF.hpp"
#include "SimdOverlapVrrRecGG.hpp"
#include "SimdOverlapVrrRecGH.hpp"
#include "SimdOverlapVrrRecGI.hpp"
#include "SimdOverlapVrrRecGP.hpp"
#include "SimdOverlapVrrRecGS.hpp"
#include "SimdOverlapVrrRecHD.hpp"
#include "SimdOverlapVrrRecHF.hpp"
#include "SimdOverlapVrrRecHG.hpp"
#include "SimdOverlapVrrRecHH.hpp"
#include "SimdOverlapVrrRecHI.hpp"
#include "SimdOverlapVrrRecHP.hpp"
#include "SimdOverlapVrrRecHS.hpp"
#include "SimdOverlapVrrRecPD.hpp"
#include "SimdOverlapVrrRecPF.hpp"
#include "SimdOverlapVrrRecPG.hpp"
#include "SimdOverlapVrrRecPH.hpp"
#include "SimdOverlapVrrRecPI.hpp"
#include "SimdOverlapVrrRecPP.hpp"
#include "SimdOverlapVrrRecPS.hpp"
#include "SimdOverlapVrrRecSD.hpp"
#include "SimdOverlapVrrRecSF.hpp"
#include "SimdOverlapVrrRecSG.hpp"
#include "SimdOverlapVrrRecSH.hpp"
#include "SimdOverlapVrrRecSI.hpp"
#include "SimdOverlapVrrRecSP.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransformHI.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_hi_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hi_kinetic_energy: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number, as every integral is a sum over them
    // and the error of a sum is bounded by the number of its terms.

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_kinetic_energy_primitive_bound, threshold / static_cast<double>(nprims));

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 4970);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 143 * nvalues, 0.0);

        return;
    }

    const auto nmax = buffer.number_of_columns();

    errors::assertMsgCritical(dimensions.size() == nprim_a * nprim_b,
                              std::string("Dimensions do not match the pairs of primitives"));

    for (size_t i = 0; i < nprim_a; i++)
    {
        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = dimensions[i * nprim_b + j];

            if (ncols == 0) continue;

            const auto p = a_exps[i] + b_exps[j];

            const auto mu = a_exps[i] * b_exps[j] / p;

            const auto fpi = mathconst::pi_value() / p;

            const auto fovl = a_norms[i] * b_norms[j] * fpi * std::sqrt(fpi);

            const auto alpha = a_exps[i];

            const auto beta = b_exps[j];

            const auto fb = a_exps[i] / p;

            const auto fa = -b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            simdfunc::compute_pb(buffer, coordinates, 3, ncols, fb);

            simdovl::compute_prim_ss_overlap(buffer, coordinates, 6, ncols, fovl, mu);

            simdovl::compute_prim_sp_overlap_0(buffer, 7, 3, 6, ncols);

            simdovl::compute_prim_sd_overlap_1(buffer, 10, 3, 6, 7, ncols, p);

            simdovl::compute_prim_sf_overlap_8(buffer, 13, 3, 7, 10, ncols, p);

            simdovl::compute_prim_sg_overlap_11(buffer, 18, 3, 10, 13, ncols, p);

            simdovl::compute_prim_sh_overlap_13(buffer, 24, 3, 13, 18, ncols, p);

            simdovl::compute_prim_si_overlap_7(buffer, 31, 3, 18, 24, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 36, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 37, 3, 7, 36, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_1(buffer, 40, 3, 6, 10, 36, 37, ncols, alpha, beta, p);

            compute_prim_sf_kinetic_energy_6(buffer, 43, 3, 7, 13, 37, 40, ncols, alpha, beta, p);

            compute_prim_sg_kinetic_energy_8(buffer, 48, 3, 10, 18, 40, 43, ncols, alpha, beta, p);

            compute_prim_sh_kinetic_energy_10(buffer, 54, 3, 13, 24, 43, 48, ncols, alpha, beta, p);

            compute_prim_si_kinetic_energy_4(buffer, 61, 3, 18, 31, 48, 54, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 66, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 69, 3, 6, 66, ncols, p);

            simdovl::compute_prim_pd_overlap_7(buffer, 72, 0, 7, 10, ncols, p);

            simdovl::compute_prim_pf_overlap_15(buffer, 75, 0, 3, 10, 13, 72, ncols, p);

            simdovl::compute_prim_pg_overlap_13(buffer, 79, 0, 3, 13, 18, 72, 75, ncols, p);

            simdovl::compute_prim_ph_overlap_9(buffer, 84, 0, 3, 18, 24, 75, 79, ncols, p);

            simdovl::compute_prim_pi_overlap_4(buffer, 90, 0, 3, 24, 31, 79, 84, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 95, 0, 36, 66, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 98, 3, 36, 69, 95, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_5(buffer, 101, 0, 37, 40, 72, ncols, alpha, beta, p);

            compute_prim_pf_kinetic_energy_16(buffer, 104, 0, 40, 43, 75, ncols, alpha, beta, p);

            compute_prim_pg_kinetic_energy_13(buffer, 107, 0, 3, 43, 48, 79, 104, ncols, alpha, beta, p);

            compute_prim_ph_kinetic_energy_9(buffer, 112, 0, 3, 48, 54, 84, 107, ncols, alpha, beta, p);

            compute_prim_pi_kinetic_energy_4(buffer, 118, 0, 54, 61, 90, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 123, 0, 6, 66, ncols, p);

            simdovl::compute_prim_dp_overlap_10(buffer, 126, 3, 66, 123, ncols, p);

            simdovl::compute_prim_dd_overlap_15(buffer, 134, 0, 3, 69, 72, 123, 126, ncols, p);

            simdovl::compute_prim_df_overlap_15(buffer, 145, 0, 3, 72, 75, 126, 134, ncols, p);

            simdovl::compute_prim_dg_overlap_11(buffer, 164, 0, 3, 75, 79, 134, 145, ncols, p);

            simdovl::compute_prim_dh_overlap_7(buffer, 192, 0, 3, 79, 84, 145, 164, ncols, p);

            simdovl::compute_prim_di_overlap_3(buffer, 220, 0, 3, 84, 90, 164, 192, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 239, 0, 6, 36, 95, 123, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_11(buffer, 242, 3, 95, 126, 239, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_15(buffer, 250, 0, 3, 98, 101, 123, 134, 239, 242, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_15(buffer, 261, 0, 3, 101, 104, 126, 145, 242, 250, ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_11(buffer, 280, 0, 3, 104, 107, 134, 164, 250, 261, ncols, alpha, beta, p);

            compute_prim_dh_kinetic_energy_7(buffer, 308, 0, 3, 107, 112, 145, 192, 261, 280, ncols, alpha, beta, p);

            compute_prim_di_kinetic_energy_3(buffer, 336, 0, 3, 112, 118, 164, 220, 280, 308, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_2(buffer, 355, 0, 66, 123, ncols, p);

            simdovl::compute_prim_fp_overlap_11(buffer, 363, 3, 123, 355, ncols, p);

            simdovl::compute_prim_fd_overlap_11(buffer, 372, 0, 3, 72, 126, 134, 355, 363, ncols, p);

            simdovl::compute_prim_ff_overlap_11(buffer, 392, 0, 3, 75, 134, 145, 363, 372, ncols, p);

            simdovl::compute_prim_fg_overlap_8(buffer, 426, 0, 3, 79, 145, 164, 372, 392, ncols, p);

            simdovl::compute_prim_fh_overlap_5(buffer, 476, 0, 3, 84, 164, 192, 392, 426, ncols, p);

            simdovl::compute_prim_fi_overlap_2(buffer, 555, 0, 3, 90, 192, 220, 426, 476, ncols, p);

            compute_prim_fs_kinetic_energy_3(buffer, 607, 0, 66, 95, 239, 355, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_11(buffer, 615, 3, 239, 363, 607, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_13(buffer, 624, 0, 3, 72, 101, 242, 250, 355, 372, 607, 615, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_11(buffer, 643, 0, 3, 75, 104, 250, 261, 363, 392, 615, 624, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_8(buffer, 676, 0, 3, 79, 107, 261, 280, 372, 426, 624, 643, ncols, alpha, beta, p);

            compute_prim_fh_kinetic_energy_5(buffer, 726, 0, 3, 84, 112, 280, 308, 392, 476, 643, 676, ncols, alpha, beta, p);

            compute_prim_fi_kinetic_energy_2(buffer, 805, 0, 3, 90, 118, 308, 336, 426, 555, 676, 726, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_7(buffer, 857, 0, 123, 355, ncols, p);

            simdovl::compute_prim_gp_overlap_7(buffer, 868, 0, 3, 355, 363, 857, ncols, p);

            simdovl::compute_prim_gd_overlap_7(buffer, 885, 0, 3, 134, 363, 372, 857, 868, ncols, p);

            simdovl::compute_prim_gf_overlap_7(buffer, 920, 0, 3, 145, 372, 392, 868, 885, ncols, p);

            simdovl::compute_prim_gg_overlap_5(buffer, 982, 0, 3, 164, 392, 426, 885, 920, ncols, p);

            simdovl::compute_prim_gh_overlap_3(buffer, 1080, 0, 3, 192, 426, 476, 920, 982, ncols, p);

            simdovl::compute_prim_gi_overlap_1(buffer, 1265, 0, 3, 220, 476, 555, 982, 1080, ncols, p);

            compute_prim_gs_kinetic_energy_4(buffer, 1419, 0, 123, 239, 607, 857, ncols, alpha, beta, p);

            compute_prim_gp_kinetic_energy_7(buffer, 1430, 0, 3, 607, 615, 868, 1419, ncols, alpha, beta, p);

            compute_prim_gd_kinetic_energy_9(buffer, 1446, 0, 3, 134, 250, 615, 624, 857, 885, 1419, 1430, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_7(buffer, 1478, 0, 3, 145, 261, 624, 643, 868, 920, 1430, 1446, ncols, alpha, beta, p);

            compute_prim_gg_kinetic_energy_5(buffer, 1536, 0, 3, 164, 280, 643, 676, 885, 982, 1446, 1478, ncols, alpha, beta, p);

            compute_prim_gh_kinetic_energy_3(buffer, 1632, 0, 3, 192, 308, 676, 726, 920, 1080, 1478, 1536, ncols, alpha, beta, p);

            compute_prim_gi_kinetic_energy_1(buffer, 1817, 0, 3, 220, 336, 726, 805, 982, 1265, 1536, 1632, ncols, alpha, beta, p);

            simdovl::compute_prim_hs_overlap_7(buffer, 1971, 0, 355, 857, ncols, p);

            simdovl::compute_prim_hp_overlap_3(buffer, 1983, 0, 3, 857, 868, 1971, ncols, p);

            simdovl::compute_prim_hd_overlap_3(buffer, 2008, 0, 3, 372, 868, 885, 1971, 1983, ncols, p);

            simdovl::compute_prim_hf_overlap_3(buffer, 2060, 0, 3, 392, 885, 920, 1983, 2008, ncols, p);

            simdovl::compute_prim_hg_overlap_2(buffer, 2156, 0, 3, 426, 920, 982, 2008, 2060, ncols, p);

            simdovl::compute_prim_hh_overlap_1(buffer, 2318, 0, 3, 476, 982, 1080, 2060, 2156, ncols, p);

            simdovl::compute_prim_hi_overlap_0(buffer, 2631, 0, 3, 555, 1080, 1265, 2156, 2318, ncols, p);

            compute_prim_hs_kinetic_energy_3(buffer, 3219, 0, 355, 607, 1419, 1971, ncols, alpha, beta, p);

            compute_prim_hp_kinetic_energy_3(buffer, 3230, 3, 1419, 1983, 3219, ncols, alpha, beta, p);

            compute_prim_hd_kinetic_energy_4(buffer, 3253, 0, 3, 372, 624, 1430, 1446, 1971, 2008, 3219, 3230, ncols, alpha, beta, p);

            compute_prim_hf_kinetic_energy_3(buffer, 3296, 0, 3, 392, 643, 1446, 1478, 1983, 2060, 3230, 3253, ncols, alpha, beta, p);

            compute_prim_hg_kinetic_energy_2(buffer, 3374, 0, 3, 426, 676, 1478, 1536, 2008, 2156, 3253, 3296, ncols, alpha, beta, p);

            compute_prim_hh_kinetic_energy_1(buffer, 3509, 0, 3, 476, 726, 1536, 1632, 2060, 2318, 3296, 3374, ncols, alpha, beta, p);

            compute_prim_hi_kinetic_energy_0(buffer, 3794, 0, 3, 555, 805, 1632, 1817, 2156, 2631, 3374, 3509, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 4382, 3794, 588, ncols);
        }
    }

    simdtrf::transform_hi(values, nvalues, buffer, 4382, nmax);

    for (size_t m = 0; m < 143; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin
