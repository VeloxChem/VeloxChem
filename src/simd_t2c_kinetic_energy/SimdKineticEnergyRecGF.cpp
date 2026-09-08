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


#include "SimdKineticEnergyRecGF.hpp"

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
#include "SimdOverlapVrrRecPD.hpp"
#include "SimdOverlapVrrRecPF.hpp"
#include "SimdOverlapVrrRecPP.hpp"
#include "SimdOverlapVrrRecPS.hpp"
#include "SimdOverlapVrrRecSD.hpp"
#include "SimdOverlapVrrRecSF.hpp"
#include "SimdOverlapVrrRecSP.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransformGF.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_gf_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gf_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 949);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 63 * nvalues, 0.0);

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

            simdovl::compute_prim_dp_overlap_3(buffer, 49, 0, 3, 24, 27, 46, ncols, p);

            simdovl::compute_prim_dd_overlap_9(buffer, 59, 0, 3, 27, 30, 46, 49, ncols, p);

            simdovl::compute_prim_df_overlap_8(buffer, 75, 0, 3, 30, 33, 49, 59, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 84, 0, 6, 15, 35, 46, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_7(buffer, 87, 0, 3, 35, 38, 49, 84, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_9(buffer, 97, 0, 3, 38, 41, 46, 59, 84, 87, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_8(buffer, 113, 0, 3, 41, 44, 49, 75, 87, 97, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_4(buffer, 122, 0, 24, 46, ncols, p);

            simdovl::compute_prim_fp_overlap_6(buffer, 131, 0, 3, 46, 49, 122, ncols, p);

            simdovl::compute_prim_fd_overlap_5(buffer, 146, 0, 3, 30, 49, 59, 122, 131, ncols, p);

            simdovl::compute_prim_ff_overlap_4(buffer, 176, 0, 3, 33, 59, 75, 131, 146, ncols, p);

            compute_prim_fs_kinetic_energy_0(buffer, 216, 0, 24, 35, 84, 122, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_6(buffer, 225, 0, 3, 84, 87, 131, 216, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_5(buffer, 240, 0, 3, 30, 41, 87, 97, 122, 146, 216, 225, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_4(buffer, 270, 0, 3, 33, 44, 97, 113, 131, 176, 225, 240, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_5(buffer, 310, 0, 46, 122, ncols, p);

            simdovl::compute_prim_gp_overlap_2(buffer, 320, 0, 3, 122, 131, 310, ncols, p);

            simdovl::compute_prim_gd_overlap_1(buffer, 344, 0, 3, 59, 131, 146, 310, 320, ncols, p);

            simdovl::compute_prim_gf_overlap_0(buffer, 406, 0, 3, 75, 146, 176, 320, 344, ncols, p);

            compute_prim_gs_kinetic_energy_1(buffer, 556, 0, 46, 84, 216, 310, ncols, alpha, beta, p);

            compute_prim_gp_kinetic_energy_2(buffer, 565, 0, 3, 216, 225, 320, 556, ncols, alpha, beta, p);

            compute_prim_gd_kinetic_energy_1(buffer, 588, 0, 3, 59, 97, 225, 240, 310, 344, 556, 565, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_0(buffer, 649, 0, 3, 75, 113, 240, 270, 320, 406, 565, 588, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 799, 649, 150, ncols);
        }
    }

    simdtrf::transform_gf(values, nvalues, buffer, 799, nmax);

    for (size_t m = 0; m < 63; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin
