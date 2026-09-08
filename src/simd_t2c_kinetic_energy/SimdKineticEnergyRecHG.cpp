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


#include "SimdKineticEnergyRecHG.hpp"

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
#include "SimdKineticEnergyVrrRecDP.hpp"
#include "SimdKineticEnergyVrrRecDS.hpp"
#include "SimdKineticEnergyVrrRecFD.hpp"
#include "SimdKineticEnergyVrrRecFF.hpp"
#include "SimdKineticEnergyVrrRecFG.hpp"
#include "SimdKineticEnergyVrrRecFP.hpp"
#include "SimdKineticEnergyVrrRecFS.hpp"
#include "SimdKineticEnergyVrrRecGD.hpp"
#include "SimdKineticEnergyVrrRecGF.hpp"
#include "SimdKineticEnergyVrrRecGG.hpp"
#include "SimdKineticEnergyVrrRecGP.hpp"
#include "SimdKineticEnergyVrrRecGS.hpp"
#include "SimdKineticEnergyVrrRecHD.hpp"
#include "SimdKineticEnergyVrrRecHF.hpp"
#include "SimdKineticEnergyVrrRecHG.hpp"
#include "SimdKineticEnergyVrrRecHP.hpp"
#include "SimdKineticEnergyVrrRecHS.hpp"
#include "SimdKineticEnergyVrrRecPD.hpp"
#include "SimdKineticEnergyVrrRecPF.hpp"
#include "SimdKineticEnergyVrrRecPG.hpp"
#include "SimdKineticEnergyVrrRecPP.hpp"
#include "SimdKineticEnergyVrrRecPS.hpp"
#include "SimdKineticEnergyVrrRecSD.hpp"
#include "SimdKineticEnergyVrrRecSF.hpp"
#include "SimdKineticEnergyVrrRecSG.hpp"
#include "SimdKineticEnergyVrrRecSP.hpp"
#include "SimdKineticEnergyVrrRecSS.hpp"
#include "SimdOverlapVrrRecDD.hpp"
#include "SimdOverlapVrrRecDF.hpp"
#include "SimdOverlapVrrRecDG.hpp"
#include "SimdOverlapVrrRecDP.hpp"
#include "SimdOverlapVrrRecDS.hpp"
#include "SimdOverlapVrrRecFD.hpp"
#include "SimdOverlapVrrRecFF.hpp"
#include "SimdOverlapVrrRecFG.hpp"
#include "SimdOverlapVrrRecFP.hpp"
#include "SimdOverlapVrrRecFS.hpp"
#include "SimdOverlapVrrRecGD.hpp"
#include "SimdOverlapVrrRecGF.hpp"
#include "SimdOverlapVrrRecGG.hpp"
#include "SimdOverlapVrrRecGP.hpp"
#include "SimdOverlapVrrRecGS.hpp"
#include "SimdOverlapVrrRecHD.hpp"
#include "SimdOverlapVrrRecHF.hpp"
#include "SimdOverlapVrrRecHG.hpp"
#include "SimdOverlapVrrRecHP.hpp"
#include "SimdOverlapVrrRecHS.hpp"
#include "SimdOverlapVrrRecPD.hpp"
#include "SimdOverlapVrrRecPF.hpp"
#include "SimdOverlapVrrRecPG.hpp"
#include "SimdOverlapVrrRecPP.hpp"
#include "SimdOverlapVrrRecPS.hpp"
#include "SimdOverlapVrrRecSD.hpp"
#include "SimdOverlapVrrRecSF.hpp"
#include "SimdOverlapVrrRecSG.hpp"
#include "SimdOverlapVrrRecSP.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransformHG.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_hg_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hg_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 2284);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 99 * nvalues, 0.0);

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

            simdovl::compute_prim_sg_overlap_9(buffer, 18, 3, 10, 13, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 21, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 22, 3, 7, 21, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_1(buffer, 25, 3, 6, 10, 21, 22, ncols, alpha, beta, p);

            compute_prim_sf_kinetic_energy_6(buffer, 28, 3, 7, 13, 22, 25, ncols, alpha, beta, p);

            compute_prim_sg_kinetic_energy_6(buffer, 33, 3, 10, 18, 25, 28, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 36, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 39, 3, 6, 36, ncols, p);

            simdovl::compute_prim_pd_overlap_7(buffer, 42, 0, 7, 10, ncols, p);

            simdovl::compute_prim_pf_overlap_15(buffer, 45, 0, 3, 10, 13, 42, ncols, p);

            simdovl::compute_prim_pg_overlap_12(buffer, 49, 0, 3, 13, 18, 42, 45, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 52, 0, 21, 36, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 55, 3, 21, 39, 52, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_5(buffer, 58, 0, 22, 25, 42, ncols, alpha, beta, p);

            compute_prim_pf_kinetic_energy_15(buffer, 61, 0, 3, 25, 28, 45, 58, ncols, alpha, beta, p);

            compute_prim_pg_kinetic_energy_9(buffer, 65, 0, 28, 33, 49, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 68, 0, 6, 36, ncols, p);

            simdovl::compute_prim_dp_overlap_10(buffer, 71, 3, 36, 68, ncols, p);

            simdovl::compute_prim_dd_overlap_14(buffer, 79, 0, 3, 39, 42, 68, 71, ncols, p);

            simdovl::compute_prim_df_overlap_13(buffer, 92, 0, 3, 42, 45, 71, 79, ncols, p);

            simdovl::compute_prim_dg_overlap_9(buffer, 104, 0, 3, 45, 49, 79, 92, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 115, 0, 6, 21, 52, 68, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_11(buffer, 118, 3, 52, 71, 115, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_14(buffer, 126, 0, 3, 55, 58, 68, 79, 115, 118, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_13(buffer, 139, 0, 3, 58, 61, 71, 92, 118, 126, ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_9(buffer, 151, 0, 3, 61, 65, 79, 104, 126, 139, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_2(buffer, 162, 0, 36, 68, ncols, p);

            simdovl::compute_prim_fp_overlap_11(buffer, 170, 3, 68, 162, ncols, p);

            simdovl::compute_prim_fd_overlap_10(buffer, 179, 0, 3, 42, 71, 79, 162, 170, ncols, p);

            simdovl::compute_prim_ff_overlap_9(buffer, 202, 0, 3, 45, 79, 92, 170, 179, ncols, p);

            simdovl::compute_prim_fg_overlap_6(buffer, 241, 0, 3, 49, 92, 104, 179, 202, ncols, p);

            compute_prim_fs_kinetic_energy_3(buffer, 268, 0, 36, 52, 115, 162, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_11(buffer, 276, 3, 115, 170, 268, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_11(buffer, 285, 0, 3, 42, 58, 118, 126, 162, 179, 268, 276, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_9(buffer, 308, 0, 3, 45, 61, 126, 139, 170, 202, 276, 285, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_6(buffer, 347, 0, 3, 49, 65, 139, 151, 179, 241, 285, 308, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_7(buffer, 374, 0, 68, 162, ncols, p);

            simdovl::compute_prim_gp_overlap_7(buffer, 385, 0, 3, 162, 170, 374, ncols, p);

            simdovl::compute_prim_gd_overlap_6(buffer, 402, 0, 3, 79, 170, 179, 374, 385, ncols, p);

            simdovl::compute_prim_gf_overlap_5(buffer, 443, 0, 3, 92, 179, 202, 385, 402, ncols, p);

            simdovl::compute_prim_gg_overlap_3(buffer, 528, 0, 3, 104, 202, 241, 402, 443, ncols, p);

            compute_prim_gs_kinetic_energy_4(buffer, 611, 0, 68, 115, 268, 374, ncols, alpha, beta, p);

            compute_prim_gp_kinetic_energy_7(buffer, 622, 0, 3, 268, 276, 385, 611, ncols, alpha, beta, p);

            compute_prim_gd_kinetic_energy_7(buffer, 638, 0, 3, 79, 126, 276, 285, 374, 402, 611, 622, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_5(buffer, 679, 0, 3, 92, 139, 285, 308, 385, 443, 622, 638, ncols, alpha, beta, p);

            compute_prim_gg_kinetic_energy_3(buffer, 764, 0, 3, 104, 151, 308, 347, 402, 528, 638, 679, ncols, alpha, beta, p);

            simdovl::compute_prim_hs_overlap_7(buffer, 847, 0, 162, 374, ncols, p);

            simdovl::compute_prim_hp_overlap_3(buffer, 859, 0, 3, 374, 385, 847, ncols, p);

            simdovl::compute_prim_hd_overlap_2(buffer, 884, 0, 3, 179, 385, 402, 847, 859, ncols, p);

            simdovl::compute_prim_hf_overlap_1(buffer, 949, 0, 3, 202, 402, 443, 859, 884, ncols, p);

            simdovl::compute_prim_hg_overlap_0(buffer, 1099, 0, 3, 241, 443, 528, 884, 949, ncols, p);

            compute_prim_hs_kinetic_energy_3(buffer, 1414, 0, 162, 268, 611, 847, ncols, alpha, beta, p);

            compute_prim_hp_kinetic_energy_3(buffer, 1425, 3, 611, 859, 1414, ncols, alpha, beta, p);

            compute_prim_hd_kinetic_energy_2(buffer, 1448, 0, 3, 179, 285, 622, 638, 847, 884, 1414, 1425, ncols, alpha, beta, p);

            compute_prim_hf_kinetic_energy_1(buffer, 1510, 0, 3, 202, 308, 638, 679, 859, 949, 1425, 1448, ncols, alpha, beta, p);

            compute_prim_hg_kinetic_energy_0(buffer, 1654, 0, 3, 241, 347, 679, 764, 884, 1099, 1448, 1510, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1969, 1654, 315, ncols);
        }
    }

    simdtrf::transform_hg(values, nvalues, buffer, 1969, nmax);

    for (size_t m = 0; m < 99; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin
