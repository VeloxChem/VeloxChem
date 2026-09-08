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
#include "SimdTransformD.hpp"
#include "SimdTransformH.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_hd_kinetic_energy(double               *values,
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 1357, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 55 * nvalues, 0.0);

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

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 16, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 17, 3, 7, 16, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_0(buffer, 20, 3, 6, 10, 16, 17, ncols, alpha, beta,
                                             p);

            simdovl::compute_prim_ps_overlap_0(buffer, 26, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_0(buffer, 29, 3, 6, 26, ncols, p);

            simdovl::compute_prim_pd_overlap_0(buffer, 38, 0, 3, 7, 10, 29, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 56, 0, 16, 26, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_0(buffer, 59, 3, 16, 29, 56, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_0(buffer, 68, 0, 3, 17, 20, 38, 59, ncols, alpha,
                                             beta, p);

            simdovl::compute_prim_ds_overlap_0(buffer, 86, 0, 6, 26, ncols, p);

            simdovl::compute_prim_dp_overlap_0(buffer, 92, 0, 3, 26, 29, 86, ncols, p);

            simdovl::compute_prim_dd_overlap_0(buffer, 110, 0, 3, 29, 38, 86, 92, ncols, p);

            compute_prim_ds_kinetic_energy_0(buffer, 146, 0, 6, 16, 56, 86, ncols, alpha, beta,
                                             p);

            compute_prim_dp_kinetic_energy_0(buffer, 152, 0, 3, 56, 59, 92, 146, ncols, alpha,
                                             beta, p);

            compute_prim_dd_kinetic_energy_0(buffer, 170, 0, 3, 59, 68, 86, 110, 146, 152, ncols,
                                             alpha, beta, p);

            simdovl::compute_prim_fs_overlap_0(buffer, 206, 0, 26, 86, ncols, p);

            simdovl::compute_prim_fp_overlap_0(buffer, 216, 0, 3, 86, 92, 206, ncols, p);

            simdovl::compute_prim_fd_overlap_0(buffer, 246, 0, 3, 38, 92, 110, 206, 216, ncols,
                                               p);

            compute_prim_fs_kinetic_energy_0(buffer, 306, 0, 26, 56, 146, 206, ncols, alpha,
                                             beta, p);

            compute_prim_fp_kinetic_energy_0(buffer, 316, 0, 3, 146, 152, 216, 306, ncols, alpha,
                                             beta, p);

            compute_prim_fd_kinetic_energy_0(buffer, 346, 0, 3, 38, 68, 152, 170, 206, 246, 306,
                                             316, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_0(buffer, 406, 0, 86, 206, ncols, p);

            simdovl::compute_prim_gp_overlap_0(buffer, 421, 0, 3, 206, 216, 406, ncols, p);

            simdovl::compute_prim_gd_overlap_0(buffer, 466, 0, 3, 110, 216, 246, 406, 421, ncols,
                                               p);

            compute_prim_gs_kinetic_energy_0(buffer, 556, 0, 86, 146, 306, 406, ncols, alpha,
                                             beta, p);

            compute_prim_gp_kinetic_energy_0(buffer, 571, 0, 3, 306, 316, 421, 556, ncols, alpha,
                                             beta, p);

            compute_prim_gd_kinetic_energy_0(buffer, 616, 0, 3, 110, 170, 316, 346, 406, 466,
                                             556, 571, ncols, alpha, beta, p);

            simdovl::compute_prim_hs_overlap_0(buffer, 706, 0, 206, 406, ncols, p);

            simdovl::compute_prim_hp_overlap_0(buffer, 727, 0, 3, 406, 421, 706, ncols, p);

            simdovl::compute_prim_hd_overlap_0(buffer, 790, 0, 3, 246, 421, 466, 706, 727, ncols,
                                               p);

            compute_prim_hs_kinetic_energy_0(buffer, 916, 0, 206, 306, 556, 706, ncols, alpha,
                                             beta, p);

            compute_prim_hp_kinetic_energy_0(buffer, 937, 0, 3, 556, 571, 727, 916, ncols, alpha,
                                             beta, p);

            compute_prim_hd_kinetic_energy_0(buffer, 1000, 0, 3, 246, 346, 571, 616, 706, 790,
                                             916, 937, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1126, 1000, 126, ncols);
        }
    }

    simdtrf::transform_d_inner(buffer, 1252, 1126, 21, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 1252, 5, nmax);

    for (size_t m = 0; m < 55; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin
