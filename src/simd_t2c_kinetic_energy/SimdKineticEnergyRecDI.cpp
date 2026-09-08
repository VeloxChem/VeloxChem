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


#include "SimdKineticEnergyRecDI.hpp"

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
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_di_kinetic_energy(double               *values,
                          const size_t          nvalues,
                          const CBasisFunction &bra,
                          const CBasisFunction &ket,
                          const CSimdMatrix    &coordinates,
                          const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_di_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 1932);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 65 * nvalues, 0.0);

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

            simdfunc::contract_primitives(buffer, 1686, 1518, 168, ncols);
        }
    }

    simdtrf::transform_i_inner(buffer, 1854, 1686, 6, nmax);

    simdtrf::transform_d_outer(values, nvalues, buffer, 1854, 13, nmax);

    for (size_t m = 0; m < 65; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin
