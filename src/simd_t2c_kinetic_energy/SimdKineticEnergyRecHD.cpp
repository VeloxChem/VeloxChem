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


#include "SimdKineticEnergyRecHD.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"

#include "SimdKineticEnergyVrrRecDD.hpp"
#include "SimdKineticEnergyVrrRecDP.hpp"
#include "SimdKineticEnergyVrrRecDS.hpp"
#include "SimdKineticEnergyVrrRecFD.hpp"
#include "SimdKineticEnergyVrrRecFP.hpp"
#include "SimdKineticEnergyVrrRecFS.hpp"
#include "SimdKineticEnergyVrrRecGD.hpp"
#include "SimdKineticEnergyVrrRecGP.hpp"
#include "SimdKineticEnergyVrrRecGS.hpp"
#include "SimdKineticEnergyVrrRecHD.hpp"
#include "SimdKineticEnergyVrrRecHP.hpp"
#include "SimdKineticEnergyVrrRecHS.hpp"
#include "SimdKineticEnergyVrrRecPD.hpp"
#include "SimdKineticEnergyVrrRecPP.hpp"
#include "SimdKineticEnergyVrrRecPS.hpp"
#include "SimdKineticEnergyVrrRecSD.hpp"
#include "SimdKineticEnergyVrrRecSP.hpp"
#include "SimdKineticEnergyVrrRecSS.hpp"
#include "SimdOverlapVrrRecDD.hpp"
#include "SimdOverlapVrrRecDP.hpp"
#include "SimdOverlapVrrRecDS.hpp"
#include "SimdOverlapVrrRecFD.hpp"
#include "SimdOverlapVrrRecFP.hpp"
#include "SimdOverlapVrrRecFS.hpp"
#include "SimdOverlapVrrRecGD.hpp"
#include "SimdOverlapVrrRecGP.hpp"
#include "SimdOverlapVrrRecGS.hpp"
#include "SimdOverlapVrrRecHD.hpp"
#include "SimdOverlapVrrRecHP.hpp"
#include "SimdOverlapVrrRecHS.hpp"
#include "SimdOverlapVrrRecPD.hpp"
#include "SimdOverlapVrrRecPP.hpp"
#include "SimdOverlapVrrRecPS.hpp"
#include "SimdOverlapVrrRecSD.hpp"
#include "SimdOverlapVrrRecSP.hpp"
#include "SimdOverlapVrrRecSS.hpp"
#include "SimdTransformHD.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_hd_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hd_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 792);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 55 * nvalues, 0.0);

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

            simdovl::compute_prim_sp_overlap_1(buffer, 7, 3, 6, ncols);

            simdovl::compute_prim_sd_overlap_3(buffer, 9, 3, 6, 7, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 11, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_1(buffer, 12, 3, 7, 11, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_2(buffer, 14, 3, 6, 9, 11, 12, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 16, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 19, 3, 6, 16, ncols, p);

            simdovl::compute_prim_pd_overlap_10(buffer, 22, 0, 7, 9, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 24, 0, 11, 16, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 27, 3, 11, 19, 24, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_8(buffer, 30, 0, 12, 14, 22, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 32, 0, 6, 16, ncols, p);

            simdovl::compute_prim_dp_overlap_8(buffer, 35, 3, 16, 32, ncols, p);

            simdovl::compute_prim_dd_overlap_12(buffer, 39, 0, 3, 19, 22, 32, 35, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 46, 0, 6, 11, 24, 32, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_9(buffer, 49, 3, 24, 35, 46, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_12(buffer, 53, 0, 3, 27, 30, 32, 39, 46, 49, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_4(buffer, 60, 0, 16, 32, ncols, p);

            simdovl::compute_prim_fp_overlap_9(buffer, 69, 0, 3, 32, 35, 60, ncols, p);

            simdovl::compute_prim_fd_overlap_8(buffer, 82, 0, 3, 22, 35, 39, 60, 69, ncols, p);

            compute_prim_fs_kinetic_energy_0(buffer, 98, 0, 16, 24, 46, 60, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_9(buffer, 107, 0, 3, 46, 49, 69, 98, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_9(buffer, 120, 0, 3, 22, 30, 49, 53, 60, 82, 98, 107, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_4(buffer, 136, 0, 32, 60, ncols, p);

            simdovl::compute_prim_gp_overlap_5(buffer, 148, 0, 3, 60, 69, 136, ncols, p);

            simdovl::compute_prim_gd_overlap_4(buffer, 172, 0, 3, 39, 69, 82, 136, 148, ncols, p);

            compute_prim_gs_kinetic_energy_0(buffer, 211, 0, 32, 46, 98, 136, ncols, alpha, beta, p);

            compute_prim_gp_kinetic_energy_5(buffer, 223, 0, 3, 98, 107, 148, 211, ncols, alpha, beta, p);

            compute_prim_gd_kinetic_energy_5(buffer, 247, 0, 3, 39, 53, 107, 120, 136, 172, 211, 223, ncols, alpha, beta, p);

            simdovl::compute_prim_hs_overlap_5(buffer, 286, 0, 60, 136, ncols, p);

            simdovl::compute_prim_hp_overlap_1(buffer, 303, 0, 3, 136, 148, 286, ncols, p);

            simdovl::compute_prim_hd_overlap_0(buffer, 350, 0, 3, 82, 148, 172, 286, 303, ncols, p);

            compute_prim_hs_kinetic_energy_1(buffer, 476, 0, 60, 98, 211, 286, ncols, alpha, beta, p);

            compute_prim_hp_kinetic_energy_1(buffer, 493, 0, 3, 211, 223, 303, 476, ncols, alpha, beta, p);

            compute_prim_hd_kinetic_energy_0(buffer, 540, 0, 3, 82, 120, 223, 247, 286, 350, 476, 493, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 666, 540, 126, ncols);
        }
    }

    simdtrf::transform_hd(values, nvalues, buffer, 666, nmax);

    for (size_t m = 0; m < 55; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin
