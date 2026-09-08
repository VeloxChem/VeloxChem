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


#include "SimdKineticEnergyRecIF.hpp"

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
#include "SimdKineticEnergyVrrRecID.hpp"
#include "SimdKineticEnergyVrrRecIF.hpp"
#include "SimdKineticEnergyVrrRecIP.hpp"
#include "SimdKineticEnergyVrrRecIS.hpp"
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
#include "SimdOverlapVrrRecID.hpp"
#include "SimdOverlapVrrRecIF.hpp"
#include "SimdOverlapVrrRecIP.hpp"
#include "SimdOverlapVrrRecIS.hpp"
#include "SimdOverlapVrrRecPD.hpp"
#include "SimdOverlapVrrRecPF.hpp"
#include "SimdOverlapVrrRecPP.hpp"
#include "SimdOverlapVrrRecPS.hpp"
#include "SimdOverlapVrrRecSD.hpp"
#include "SimdOverlapVrrRecSF.hpp"
#include "SimdOverlapVrrRecSP.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransformIF.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_if_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_if_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 1933);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 91 * nvalues, 0.0);

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

            simdovl::compute_prim_dp_overlap_8(buffer, 49, 3, 24, 46, ncols, p);

            simdovl::compute_prim_dd_overlap_16(buffer, 53, 0, 3, 27, 30, 46, 49, ncols, p);

            simdovl::compute_prim_df_overlap_16(buffer, 60, 0, 3, 30, 33, 49, 53, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 67, 0, 6, 15, 35, 46, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_9(buffer, 70, 3, 35, 49, 67, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_16(buffer, 74, 0, 3, 38, 41, 46, 53, 67, 70, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_16(buffer, 81, 0, 3, 41, 44, 49, 60, 70, 74, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_2(buffer, 88, 0, 24, 46, ncols, p);

            simdovl::compute_prim_fp_overlap_13(buffer, 96, 0, 3, 46, 49, 88, ncols, p);

            simdovl::compute_prim_fd_overlap_13(buffer, 105, 0, 3, 30, 49, 53, 88, 96, ncols, p);

            simdovl::compute_prim_ff_overlap_12(buffer, 121, 0, 3, 33, 53, 60, 96, 105, ncols, p);

            compute_prim_fs_kinetic_energy_3(buffer, 135, 0, 24, 35, 67, 88, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_13(buffer, 143, 0, 3, 67, 70, 96, 135, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_15(buffer, 152, 0, 3, 30, 41, 70, 74, 88, 105, 135, 143, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_12(buffer, 168, 0, 3, 33, 44, 74, 81, 96, 121, 143, 152, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_9(buffer, 182, 0, 46, 88, ncols, p);

            simdovl::compute_prim_gp_overlap_10(buffer, 194, 0, 3, 88, 96, 182, ncols, p);

            simdovl::compute_prim_gd_overlap_9(buffer, 210, 0, 3, 53, 96, 105, 182, 194, ncols, p);

            simdovl::compute_prim_gf_overlap_8(buffer, 243, 0, 3, 60, 105, 121, 194, 210, ncols, p);

            compute_prim_gs_kinetic_energy_7(buffer, 267, 0, 46, 67, 135, 182, ncols, alpha, beta, p);

            compute_prim_gp_kinetic_energy_10(buffer, 279, 0, 3, 135, 143, 194, 267, ncols, alpha, beta, p);

            compute_prim_gd_kinetic_energy_11(buffer, 295, 0, 3, 53, 74, 143, 152, 182, 210, 267, 279, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_8(buffer, 328, 0, 3, 60, 81, 152, 168, 194, 243, 279, 295, ncols, alpha, beta, p);

            simdovl::compute_prim_hs_overlap_9(buffer, 352, 0, 88, 182, ncols, p);

            simdovl::compute_prim_hp_overlap_6(buffer, 369, 0, 3, 182, 194, 352, ncols, p);

            simdovl::compute_prim_hd_overlap_5(buffer, 399, 0, 3, 105, 194, 210, 352, 369, ncols, p);

            simdovl::compute_prim_hf_overlap_4(buffer, 469, 0, 3, 121, 210, 243, 369, 399, ncols, p);

            compute_prim_hs_kinetic_energy_6(buffer, 542, 0, 88, 135, 267, 352, ncols, alpha, beta, p);

            compute_prim_hp_kinetic_energy_6(buffer, 559, 0, 3, 267, 279, 369, 542, ncols, alpha, beta, p);

            compute_prim_hd_kinetic_energy_6(buffer, 589, 0, 3, 105, 152, 279, 295, 352, 399, 542, 559, ncols, alpha, beta, p);

            compute_prim_hf_kinetic_energy_4(buffer, 659, 0, 3, 121, 168, 295, 328, 369, 469, 559, 589, ncols, alpha, beta, p);

            simdovl::compute_prim_is_overlap_5(buffer, 732, 0, 182, 352, ncols, p);

            simdovl::compute_prim_ip_overlap_2(buffer, 751, 0, 3, 352, 369, 732, ncols, p);

            simdovl::compute_prim_id_overlap_1(buffer, 794, 0, 3, 210, 369, 399, 732, 751, ncols, p);

            simdovl::compute_prim_if_overlap_0(buffer, 914, 0, 3, 243, 399, 469, 751, 794, ncols, p);

            compute_prim_is_kinetic_energy_2(buffer, 1194, 0, 182, 267, 542, 732, ncols, alpha, beta, p);

            compute_prim_ip_kinetic_energy_2(buffer, 1212, 0, 3, 542, 559, 751, 1194, ncols, alpha, beta, p);

            compute_prim_id_kinetic_energy_1(buffer, 1254, 0, 3, 210, 295, 559, 589, 732, 794, 1194, 1212, ncols, alpha, beta, p);

            compute_prim_if_kinetic_energy_0(buffer, 1373, 0, 3, 243, 328, 589, 659, 751, 914, 1212, 1254, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1653, 1373, 280, ncols);
        }
    }

    simdtrf::transform_if(values, nvalues, buffer, 1653, nmax);

    for (size_t m = 0; m < 91; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin
