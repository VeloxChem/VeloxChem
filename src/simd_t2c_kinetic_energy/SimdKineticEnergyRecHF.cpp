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


#include "SimdKineticEnergyRecHF.hpp"

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
#include "SimdKineticEnergyVrrRecDP.hpp"
#include "SimdKineticEnergyVrrRecDS.hpp"
#include "SimdKineticEnergyVrrRecFD.hpp"
#include "SimdKineticEnergyVrrRecFF.hpp"
#include "SimdKineticEnergyVrrRecFP.hpp"
#include "SimdKineticEnergyVrrRecFS.hpp"
#include "SimdKineticEnergyVrrRecGD.hpp"
#include "SimdKineticEnergyVrrRecGF.hpp"
#include "SimdKineticEnergyVrrRecGP.hpp"
#include "SimdKineticEnergyVrrRecGS.hpp"
#include "SimdKineticEnergyVrrRecHD.hpp"
#include "SimdKineticEnergyVrrRecHF.hpp"
#include "SimdKineticEnergyVrrRecHP.hpp"
#include "SimdKineticEnergyVrrRecHS.hpp"
#include "SimdKineticEnergyVrrRecPD.hpp"
#include "SimdKineticEnergyVrrRecPF.hpp"
#include "SimdKineticEnergyVrrRecPP.hpp"
#include "SimdKineticEnergyVrrRecPS.hpp"
#include "SimdKineticEnergyVrrRecSD.hpp"
#include "SimdKineticEnergyVrrRecSF.hpp"
#include "SimdKineticEnergyVrrRecSP.hpp"
#include "SimdKineticEnergyVrrRecSS.hpp"
#include "SimdOverlapVrrRecDD.hpp"
#include "SimdOverlapVrrRecDF.hpp"
#include "SimdOverlapVrrRecDP.hpp"
#include "SimdOverlapVrrRecDS.hpp"
#include "SimdOverlapVrrRecFD.hpp"
#include "SimdOverlapVrrRecFF.hpp"
#include "SimdOverlapVrrRecFP.hpp"
#include "SimdOverlapVrrRecFS.hpp"
#include "SimdOverlapVrrRecGD.hpp"
#include "SimdOverlapVrrRecGF.hpp"
#include "SimdOverlapVrrRecGP.hpp"
#include "SimdOverlapVrrRecGS.hpp"
#include "SimdOverlapVrrRecHD.hpp"
#include "SimdOverlapVrrRecHF.hpp"
#include "SimdOverlapVrrRecHP.hpp"
#include "SimdOverlapVrrRecHS.hpp"
#include "SimdOverlapVrrRecPD.hpp"
#include "SimdOverlapVrrRecPF.hpp"
#include "SimdOverlapVrrRecPP.hpp"
#include "SimdOverlapVrrRecPS.hpp"
#include "SimdOverlapVrrRecSD.hpp"
#include "SimdOverlapVrrRecSF.hpp"
#include "SimdOverlapVrrRecSP.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransformHF.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_hf_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hf_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 1391);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 77 * nvalues, 0.0);

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

            simdovl::compute_prim_sf_overlap_7(buffer, 13, 3, 7, 10, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 15, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 16, 3, 7, 15, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_1(buffer, 19, 3, 6, 10, 15, 16, ncols, alpha, beta, p);

            compute_prim_sf_kinetic_energy_5(buffer, 22, 3, 7, 13, 16, 19, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 24, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 27, 3, 6, 24, ncols, p);

            simdovl::compute_prim_pd_overlap_7(buffer, 30, 0, 7, 10, ncols, p);

            simdovl::compute_prim_pf_overlap_11(buffer, 33, 0, 10, 13, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 35, 0, 15, 24, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 38, 3, 15, 27, 35, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_5(buffer, 41, 0, 16, 19, 30, ncols, alpha, beta, p);

            compute_prim_pf_kinetic_energy_11(buffer, 44, 0, 19, 22, 33, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 46, 0, 6, 24, ncols, p);

            simdovl::compute_prim_dp_overlap_9(buffer, 49, 0, 3, 24, 27, 46, ncols, p);

            simdovl::compute_prim_dd_overlap_13(buffer, 56, 0, 3, 27, 30, 46, 49, ncols, p);

            simdovl::compute_prim_df_overlap_12(buffer, 65, 0, 3, 30, 33, 49, 56, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 72, 0, 6, 15, 35, 46, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_10(buffer, 75, 0, 3, 35, 38, 49, 72, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_13(buffer, 82, 0, 3, 38, 41, 46, 56, 72, 75, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_12(buffer, 91, 0, 3, 41, 44, 49, 65, 75, 82, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_4(buffer, 98, 0, 24, 46, ncols, p);

            simdovl::compute_prim_fp_overlap_10(buffer, 107, 0, 3, 46, 49, 98, ncols, p);

            simdovl::compute_prim_fd_overlap_9(buffer, 119, 0, 3, 30, 49, 56, 98, 107, ncols, p);

            simdovl::compute_prim_ff_overlap_8(buffer, 142, 0, 3, 33, 56, 65, 107, 119, ncols, p);

            compute_prim_fs_kinetic_energy_0(buffer, 158, 0, 24, 35, 72, 98, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_10(buffer, 167, 0, 3, 72, 75, 107, 158, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_10(buffer, 179, 0, 3, 30, 41, 75, 82, 98, 119, 158, 167, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_8(buffer, 202, 0, 3, 33, 44, 82, 91, 107, 142, 167, 179, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_4(buffer, 218, 0, 46, 98, ncols, p);

            simdovl::compute_prim_gp_overlap_6(buffer, 230, 0, 3, 98, 107, 218, ncols, p);

            simdovl::compute_prim_gd_overlap_5(buffer, 252, 0, 3, 56, 107, 119, 218, 230, ncols, p);

            simdovl::compute_prim_gf_overlap_4(buffer, 300, 0, 3, 65, 119, 142, 230, 252, ncols, p);

            compute_prim_gs_kinetic_energy_0(buffer, 355, 0, 46, 72, 158, 218, ncols, alpha, beta, p);

            compute_prim_gp_kinetic_energy_6(buffer, 367, 0, 3, 158, 167, 230, 355, ncols, alpha, beta, p);

            compute_prim_gd_kinetic_energy_6(buffer, 389, 0, 3, 56, 82, 167, 179, 218, 252, 355, 367, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_4(buffer, 437, 0, 3, 65, 91, 179, 202, 230, 300, 367, 389, ncols, alpha, beta, p);

            simdovl::compute_prim_hs_overlap_6(buffer, 492, 0, 98, 218, ncols, p);

            simdovl::compute_prim_hp_overlap_2(buffer, 506, 0, 3, 218, 230, 492, ncols, p);

            simdovl::compute_prim_hd_overlap_1(buffer, 539, 0, 3, 119, 230, 252, 492, 506, ncols, p);

            simdovl::compute_prim_hf_overlap_0(buffer, 628, 0, 3, 142, 252, 300, 506, 539, ncols, p);

            compute_prim_hs_kinetic_energy_2(buffer, 838, 0, 98, 158, 355, 492, ncols, alpha, beta, p);

            compute_prim_hp_kinetic_energy_2(buffer, 851, 0, 3, 355, 367, 506, 838, ncols, alpha, beta, p);

            compute_prim_hd_kinetic_energy_1(buffer, 883, 0, 3, 119, 179, 367, 389, 492, 539, 838, 851, ncols, alpha, beta, p);

            compute_prim_hf_kinetic_energy_0(buffer, 971, 0, 3, 142, 202, 389, 437, 506, 628, 851, 883, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1181, 971, 210, ncols);
        }
    }

    simdtrf::transform_hf(values, nvalues, buffer, 1181, nmax);

    for (size_t m = 0; m < 77; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin
