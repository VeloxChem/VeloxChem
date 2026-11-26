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
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_hi_kinetic_energy(double               *values,
                          const size_t          nvalues,
                          const CBasisFunction &bra,
                          const CBasisFunction &ket,
                          const CSimdMatrix    &coordinates,
                          CSimdMatrix          &buffer,
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 10275, 9414, 588, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 143 * nvalues, 0.0);

        return;
    }

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

            simdovl::compute_prim_sd_overlap_0(buffer, 10, 3, 6, 7, ncols, p);

            simdovl::compute_prim_sf_overlap_0(buffer, 16, 3, 7, 10, ncols, p);

            simdovl::compute_prim_sg_overlap_0(buffer, 26, 3, 10, 16, ncols, p);

            simdovl::compute_prim_sh_overlap_0(buffer, 41, 3, 16, 26, ncols, p);

            simdovl::compute_prim_si_overlap_0(buffer, 62, 3, 26, 41, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 90, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 91, 3, 7, 90, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_0(buffer, 94, 3, 6, 10, 90, 91, ncols, alpha, beta,
                                             p);

            compute_prim_sf_kinetic_energy_0(buffer, 100, 3, 7, 16, 91, 94, ncols, alpha, beta,
                                             p);

            compute_prim_sg_kinetic_energy_0(buffer, 110, 3, 10, 26, 94, 100, ncols, alpha, beta,
                                             p);

            compute_prim_sh_kinetic_energy_0(buffer, 125, 3, 16, 41, 100, 110, ncols, alpha,
                                             beta, p);

            compute_prim_si_kinetic_energy_0(buffer, 146, 3, 26, 62, 110, 125, ncols, alpha,
                                             beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 174, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_0(buffer, 177, 3, 6, 174, ncols, p);

            simdovl::compute_prim_pd_overlap_0(buffer, 186, 0, 3, 7, 10, 177, ncols, p);

            simdovl::compute_prim_pf_overlap_0(buffer, 204, 0, 3, 10, 16, 177, 186, ncols, p);

            simdovl::compute_prim_pg_overlap_0(buffer, 234, 0, 3, 16, 26, 186, 204, ncols, p);

            simdovl::compute_prim_ph_overlap_0(buffer, 279, 0, 3, 26, 41, 204, 234, ncols, p);

            simdovl::compute_prim_pi_overlap_0(buffer, 342, 0, 3, 41, 62, 234, 279, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 426, 0, 90, 174, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_0(buffer, 429, 3, 90, 177, 426, ncols, alpha, beta,
                                             p);

            compute_prim_pd_kinetic_energy_0(buffer, 438, 0, 3, 91, 94, 186, 429, ncols, alpha,
                                             beta, p);

            compute_prim_pf_kinetic_energy_0(buffer, 456, 0, 3, 94, 100, 204, 438, ncols, alpha,
                                             beta, p);

            compute_prim_pg_kinetic_energy_0(buffer, 486, 0, 3, 100, 110, 234, 456, ncols, alpha,
                                             beta, p);

            compute_prim_ph_kinetic_energy_0(buffer, 531, 0, 3, 110, 125, 279, 486, ncols, alpha,
                                             beta, p);

            compute_prim_pi_kinetic_energy_0(buffer, 594, 0, 3, 125, 146, 342, 531, ncols, alpha,
                                             beta, p);

            simdovl::compute_prim_ds_overlap_0(buffer, 678, 0, 6, 174, ncols, p);

            simdovl::compute_prim_dp_overlap_0(buffer, 684, 0, 3, 174, 177, 678, ncols, p);

            simdovl::compute_prim_dd_overlap_0(buffer, 702, 0, 3, 177, 186, 678, 684, ncols, p);

            simdovl::compute_prim_df_overlap_0(buffer, 738, 0, 3, 186, 204, 684, 702, ncols, p);

            simdovl::compute_prim_dg_overlap_0(buffer, 798, 0, 3, 204, 234, 702, 738, ncols, p);

            simdovl::compute_prim_dh_overlap_0(buffer, 888, 0, 3, 234, 279, 738, 798, ncols, p);

            simdovl::compute_prim_di_overlap_0(buffer, 1014, 0, 3, 279, 342, 798, 888, ncols,
                                               p);

            compute_prim_ds_kinetic_energy_0(buffer, 1182, 0, 6, 90, 426, 678, ncols, alpha,
                                             beta, p);

            compute_prim_dp_kinetic_energy_0(buffer, 1188, 0, 3, 426, 429, 684, 1182, ncols,
                                             alpha, beta, p);

            compute_prim_dd_kinetic_energy_0(buffer, 1206, 0, 3, 429, 438, 678, 702, 1182, 1188,
                                             ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_0(buffer, 1242, 0, 3, 438, 456, 684, 738, 1188, 1206,
                                             ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_0(buffer, 1302, 0, 3, 456, 486, 702, 798, 1206, 1242,
                                             ncols, alpha, beta, p);

            compute_prim_dh_kinetic_energy_0(buffer, 1392, 0, 3, 486, 531, 738, 888, 1242, 1302,
                                             ncols, alpha, beta, p);

            compute_prim_di_kinetic_energy_0(buffer, 1518, 0, 3, 531, 594, 798, 1014, 1302, 1392,
                                             ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_0(buffer, 1686, 0, 174, 678, ncols, p);

            simdovl::compute_prim_fp_overlap_0(buffer, 1696, 0, 3, 678, 684, 1686, ncols, p);

            simdovl::compute_prim_fd_overlap_0(buffer, 1726, 0, 3, 186, 684, 702, 1686, 1696,
                                               ncols, p);

            simdovl::compute_prim_ff_overlap_0(buffer, 1786, 0, 3, 204, 702, 738, 1696, 1726,
                                               ncols, p);

            simdovl::compute_prim_fg_overlap_0(buffer, 1886, 0, 3, 234, 738, 798, 1726, 1786,
                                               ncols, p);

            simdovl::compute_prim_fh_overlap_0(buffer, 2036, 0, 3, 279, 798, 888, 1786, 1886,
                                               ncols, p);

            simdovl::compute_prim_fi_overlap_0(buffer, 2246, 0, 3, 342, 888, 1014, 1886, 2036,
                                               ncols, p);

            compute_prim_fs_kinetic_energy_0(buffer, 2526, 0, 174, 426, 1182, 1686, ncols, alpha,
                                             beta, p);

            compute_prim_fp_kinetic_energy_0(buffer, 2536, 0, 3, 1182, 1188, 1696, 2526, ncols,
                                             alpha, beta, p);

            compute_prim_fd_kinetic_energy_0(buffer, 2566, 0, 3, 186, 438, 1188, 1206, 1686,
                                             1726, 2526, 2536, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_0(buffer, 2626, 0, 3, 204, 456, 1206, 1242, 1696,
                                             1786, 2536, 2566, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_0(buffer, 2726, 0, 3, 234, 486, 1242, 1302, 1726,
                                             1886, 2566, 2626, ncols, alpha, beta, p);

            compute_prim_fh_kinetic_energy_0(buffer, 2876, 0, 3, 279, 531, 1302, 1392, 1786,
                                             2036, 2626, 2726, ncols, alpha, beta, p);

            compute_prim_fi_kinetic_energy_0(buffer, 3086, 0, 3, 342, 594, 1392, 1518, 1886,
                                             2246, 2726, 2876, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_0(buffer, 3366, 0, 678, 1686, ncols, p);

            simdovl::compute_prim_gp_overlap_0(buffer, 3381, 0, 3, 1686, 1696, 3366, ncols, p);

            simdovl::compute_prim_gd_overlap_0(buffer, 3426, 0, 3, 702, 1696, 1726, 3366, 3381,
                                               ncols, p);

            simdovl::compute_prim_gf_overlap_0(buffer, 3516, 0, 3, 738, 1726, 1786, 3381, 3426,
                                               ncols, p);

            simdovl::compute_prim_gg_overlap_0(buffer, 3666, 0, 3, 798, 1786, 1886, 3426, 3516,
                                               ncols, p);

            simdovl::compute_prim_gh_overlap_0(buffer, 3891, 0, 3, 888, 1886, 2036, 3516, 3666,
                                               ncols, p);

            simdovl::compute_prim_gi_overlap_0(buffer, 4206, 0, 3, 1014, 2036, 2246, 3666, 3891,
                                               ncols, p);

            compute_prim_gs_kinetic_energy_0(buffer, 4626, 0, 678, 1182, 2526, 3366, ncols,
                                             alpha, beta, p);

            compute_prim_gp_kinetic_energy_0(buffer, 4641, 0, 3, 2526, 2536, 3381, 4626, ncols,
                                             alpha, beta, p);

            compute_prim_gd_kinetic_energy_0(buffer, 4686, 0, 3, 702, 1206, 2536, 2566, 3366,
                                             3426, 4626, 4641, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_0(buffer, 4776, 0, 3, 738, 1242, 2566, 2626, 3381,
                                             3516, 4641, 4686, ncols, alpha, beta, p);

            compute_prim_gg_kinetic_energy_0(buffer, 4926, 0, 3, 798, 1302, 2626, 2726, 3426,
                                             3666, 4686, 4776, ncols, alpha, beta, p);

            compute_prim_gh_kinetic_energy_0(buffer, 5151, 0, 3, 888, 1392, 2726, 2876, 3516,
                                             3891, 4776, 4926, ncols, alpha, beta, p);

            compute_prim_gi_kinetic_energy_0(buffer, 5466, 0, 3, 1014, 1518, 2876, 3086, 3666,
                                             4206, 4926, 5151, ncols, alpha, beta, p);

            simdovl::compute_prim_hs_overlap_0(buffer, 5886, 0, 1686, 3366, ncols, p);

            simdovl::compute_prim_hp_overlap_0(buffer, 5907, 0, 3, 3366, 3381, 5886, ncols, p);

            simdovl::compute_prim_hd_overlap_0(buffer, 5970, 0, 3, 1726, 3381, 3426, 5886, 5907,
                                               ncols, p);

            simdovl::compute_prim_hf_overlap_0(buffer, 6096, 0, 3, 1786, 3426, 3516, 5907, 5970,
                                               ncols, p);

            simdovl::compute_prim_hg_overlap_0(buffer, 6306, 0, 3, 1886, 3516, 3666, 5970, 6096,
                                               ncols, p);

            simdovl::compute_prim_hh_overlap_0(buffer, 6621, 0, 3, 2036, 3666, 3891, 6096, 6306,
                                               ncols, p);

            simdovl::compute_prim_hi_overlap_0(buffer, 7062, 0, 3, 2246, 3891, 4206, 6306, 6621,
                                               ncols, p);

            compute_prim_hs_kinetic_energy_0(buffer, 7650, 0, 1686, 2526, 4626, 5886, ncols,
                                             alpha, beta, p);

            compute_prim_hp_kinetic_energy_0(buffer, 7671, 0, 3, 4626, 4641, 5907, 7650, ncols,
                                             alpha, beta, p);

            compute_prim_hd_kinetic_energy_0(buffer, 7734, 0, 3, 1726, 2566, 4641, 4686, 5886,
                                             5970, 7650, 7671, ncols, alpha, beta, p);

            compute_prim_hf_kinetic_energy_0(buffer, 7860, 0, 3, 1786, 2626, 4686, 4776, 5907,
                                             6096, 7671, 7734, ncols, alpha, beta, p);

            compute_prim_hg_kinetic_energy_0(buffer, 8070, 0, 3, 1886, 2726, 4776, 4926, 5970,
                                             6306, 7734, 7860, ncols, alpha, beta, p);

            compute_prim_hh_kinetic_energy_0(buffer, 8385, 0, 3, 2036, 2876, 4926, 5151, 6096,
                                             6621, 7860, 8070, ncols, alpha, beta, p);

            compute_prim_hi_kinetic_energy_0(buffer, 8826, 0, 3, 2246, 3086, 5151, 5466, 6306,
                                             7062, 8070, 8385, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 9414, 8826, 588, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 10002, 9414, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 10002, 13, nmax);

    for (size_t m = 0; m < 143; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin
