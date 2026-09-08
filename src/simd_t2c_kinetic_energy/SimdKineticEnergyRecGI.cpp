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


#include "SimdKineticEnergyRecGI.hpp"

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
#include "SimdTransformGI.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_gi_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_gi_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 3410);

    if (buffer.number_of_columns() == 0)
    {
        std::fill(values, values + 117 * nvalues, 0.0);

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

            simdovl::compute_prim_sf_overlap_5(buffer, 13, 3, 7, 10, ncols, p);

            simdovl::compute_prim_sg_overlap_6(buffer, 22, 3, 10, 13, ncols, p);

            simdovl::compute_prim_sh_overlap_11(buffer, 34, 3, 13, 22, ncols, p);

            simdovl::compute_prim_si_overlap_6(buffer, 45, 3, 22, 34, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 50, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 51, 3, 7, 50, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_1(buffer, 54, 3, 6, 10, 50, 51, ncols, alpha, beta, p);

            compute_prim_sf_kinetic_energy_3(buffer, 57, 3, 7, 13, 51, 54, ncols, alpha, beta, p);

            compute_prim_sg_kinetic_energy_3(buffer, 66, 3, 10, 22, 54, 57, ncols, alpha, beta, p);

            compute_prim_sh_kinetic_energy_8(buffer, 78, 3, 13, 34, 57, 66, ncols, alpha, beta, p);

            compute_prim_si_kinetic_energy_3(buffer, 89, 3, 22, 45, 66, 78, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 94, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 97, 3, 6, 94, ncols, p);

            simdovl::compute_prim_pd_overlap_7(buffer, 100, 0, 7, 10, ncols, p);

            simdovl::compute_prim_pf_overlap_14(buffer, 103, 0, 3, 10, 13, 100, ncols, p);

            simdovl::compute_prim_pg_overlap_11(buffer, 109, 0, 3, 13, 22, 100, 103, ncols, p);

            simdovl::compute_prim_ph_overlap_7(buffer, 119, 0, 3, 22, 34, 103, 109, ncols, p);

            simdovl::compute_prim_pi_overlap_3(buffer, 129, 0, 3, 34, 45, 109, 119, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 134, 0, 50, 94, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 137, 3, 50, 97, 134, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_5(buffer, 140, 0, 51, 54, 100, ncols, alpha, beta, p);

            compute_prim_pf_kinetic_energy_14(buffer, 143, 0, 54, 57, 103, ncols, alpha, beta, p);

            compute_prim_pg_kinetic_energy_11(buffer, 146, 0, 3, 57, 66, 109, 143, ncols, alpha, beta, p);

            compute_prim_ph_kinetic_energy_7(buffer, 156, 0, 3, 66, 78, 119, 146, ncols, alpha, beta, p);

            compute_prim_pi_kinetic_energy_3(buffer, 166, 0, 78, 89, 129, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 171, 0, 6, 94, ncols, p);

            simdovl::compute_prim_dp_overlap_7(buffer, 174, 3, 94, 171, ncols, p);

            simdovl::compute_prim_dd_overlap_11(buffer, 183, 0, 3, 97, 100, 171, 174, ncols, p);

            simdovl::compute_prim_df_overlap_11(buffer, 196, 0, 3, 100, 103, 174, 183, ncols, p);

            simdovl::compute_prim_dg_overlap_8(buffer, 218, 0, 3, 103, 109, 183, 196, ncols, p);

            simdovl::compute_prim_dh_overlap_5(buffer, 250, 0, 3, 109, 119, 196, 218, ncols, p);

            simdovl::compute_prim_di_overlap_2(buffer, 306, 0, 3, 119, 129, 218, 250, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 342, 0, 6, 50, 134, 171, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_8(buffer, 345, 3, 134, 174, 342, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_11(buffer, 354, 0, 3, 137, 140, 171, 183, 342, 345, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_11(buffer, 367, 0, 3, 140, 143, 174, 196, 345, 354, ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_8(buffer, 389, 0, 3, 143, 146, 183, 218, 354, 367, ncols, alpha, beta, p);

            compute_prim_dh_kinetic_energy_5(buffer, 421, 0, 3, 146, 156, 196, 250, 367, 389, ncols, alpha, beta, p);

            compute_prim_di_kinetic_energy_2(buffer, 477, 0, 3, 156, 166, 218, 306, 389, 421, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_2(buffer, 513, 0, 94, 171, ncols, p);

            simdovl::compute_prim_fp_overlap_7(buffer, 521, 0, 3, 171, 174, 513, ncols, p);

            simdovl::compute_prim_fd_overlap_7(buffer, 533, 0, 3, 100, 174, 183, 513, 521, ncols, p);

            simdovl::compute_prim_ff_overlap_7(buffer, 558, 0, 3, 103, 183, 196, 521, 533, ncols, p);

            simdovl::compute_prim_fg_overlap_5(buffer, 601, 0, 3, 109, 196, 218, 533, 558, ncols, p);

            simdovl::compute_prim_fh_overlap_3(buffer, 666, 0, 3, 119, 218, 250, 558, 601, ncols, p);

            simdovl::compute_prim_fi_overlap_1(buffer, 787, 0, 3, 129, 250, 306, 601, 666, ncols, p);

            compute_prim_fs_kinetic_energy_3(buffer, 896, 0, 94, 134, 342, 513, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_7(buffer, 904, 0, 3, 342, 345, 521, 896, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_8(buffer, 915, 0, 3, 100, 140, 345, 354, 513, 533, 896, 904, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_7(buffer, 937, 0, 3, 103, 143, 354, 367, 521, 558, 904, 915, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_5(buffer, 976, 0, 3, 109, 146, 367, 389, 533, 601, 915, 937, ncols, alpha, beta, p);

            compute_prim_fh_kinetic_energy_3(buffer, 1039, 0, 3, 119, 156, 389, 421, 558, 666, 937, 976, ncols, alpha, beta, p);

            compute_prim_fi_kinetic_energy_1(buffer, 1160, 0, 3, 129, 166, 421, 477, 601, 787, 976, 1039, ncols, alpha, beta, p);

            simdovl::compute_prim_gs_overlap_6(buffer, 1269, 0, 171, 513, ncols, p);

            simdovl::compute_prim_gp_overlap_3(buffer, 1278, 0, 3, 513, 521, 1269, ncols, p);

            simdovl::compute_prim_gd_overlap_3(buffer, 1297, 0, 3, 183, 521, 533, 1269, 1278, ncols, p);

            simdovl::compute_prim_gf_overlap_3(buffer, 1336, 0, 3, 196, 533, 558, 1278, 1297, ncols, p);

            simdovl::compute_prim_gg_overlap_2(buffer, 1408, 0, 3, 218, 558, 601, 1297, 1336, ncols, p);

            simdovl::compute_prim_gh_overlap_1(buffer, 1528, 0, 3, 250, 601, 666, 1336, 1408, ncols, p);

            simdovl::compute_prim_gi_overlap_0(buffer, 1752, 0, 3, 306, 666, 787, 1408, 1528, ncols, p);

            compute_prim_gs_kinetic_energy_2(buffer, 2172, 0, 171, 342, 896, 1269, ncols, alpha, beta, p);

            compute_prim_gp_kinetic_energy_3(buffer, 2180, 3, 896, 1278, 2172, ncols, alpha, beta, p);

            compute_prim_gd_kinetic_energy_4(buffer, 2197, 0, 3, 183, 354, 904, 915, 1269, 1297, 2172, 2180, ncols, alpha, beta, p);

            compute_prim_gf_kinetic_energy_3(buffer, 2227, 0, 3, 196, 367, 915, 937, 1278, 1336, 2180, 2197, ncols, alpha, beta, p);

            compute_prim_gg_kinetic_energy_2(buffer, 2281, 0, 3, 218, 389, 937, 976, 1297, 1408, 2197, 2227, ncols, alpha, beta, p);

            compute_prim_gh_kinetic_energy_1(buffer, 2374, 0, 3, 250, 421, 976, 1039, 1336, 1528, 2227, 2281, ncols, alpha, beta, p);

            compute_prim_gi_kinetic_energy_0(buffer, 2570, 0, 3, 306, 477, 1039, 1160, 1408, 1752, 2281, 2374, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 2990, 2570, 420, ncols);
        }
    }

    simdtrf::transform_gi(values, nvalues, buffer, 2990, nmax);

    for (size_t m = 0; m < 117; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin
