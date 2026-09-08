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


#include "SimdKineticEnergyRecFI.hpp"

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
#include "SimdTransformFI.hpp"

namespace simdkin {  // simdkin namespace

auto
compute_fi_kinetic_energy(double               *values,
                               const size_t          nvalues,
                               const CBasisFunction &bra,
                               const CBasisFunction &ket,
                               const CSimdMatrix    &coordinates,
                               const double          threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_fi_kinetic_energy: Number of values exceeds number of atom pairs"));
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

    auto buffer = simdfunc::make_primitive_buffer(dimensions, 2144);

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

            simdovl::compute_prim_sf_overlap_5(buffer, 13, 3, 7, 10, ncols, p);

            simdovl::compute_prim_sg_overlap_6(buffer, 22, 3, 10, 13, ncols, p);

            simdovl::compute_prim_sh_overlap_9(buffer, 34, 3, 13, 22, ncols, p);

            simdovl::compute_prim_si_overlap_5(buffer, 51, 3, 22, 34, ncols, p);

            compute_prim_ss_kinetic_energy_0(buffer, coordinates, 64, 6, ncols, mu);

            compute_prim_sp_kinetic_energy_0(buffer, 65, 3, 7, 64, ncols, alpha, beta, p);

            compute_prim_sd_kinetic_energy_1(buffer, 68, 3, 6, 10, 64, 65, ncols, alpha, beta, p);

            compute_prim_sf_kinetic_energy_3(buffer, 71, 3, 7, 13, 65, 68, ncols, alpha, beta, p);

            compute_prim_sg_kinetic_energy_3(buffer, 80, 3, 10, 22, 68, 71, ncols, alpha, beta, p);

            compute_prim_sh_kinetic_energy_6(buffer, 92, 3, 13, 34, 71, 80, ncols, alpha, beta, p);

            compute_prim_si_kinetic_energy_2(buffer, 109, 3, 22, 51, 80, 92, ncols, alpha, beta, p);

            simdovl::compute_prim_ps_overlap_0(buffer, 122, 0, 6, ncols);

            simdovl::compute_prim_pp_overlap_2(buffer, 125, 3, 6, 122, ncols, p);

            simdovl::compute_prim_pd_overlap_9(buffer, 128, 0, 3, 7, 10, 125, ncols, p);

            simdovl::compute_prim_pf_overlap_10(buffer, 133, 0, 3, 10, 13, 125, 128, ncols, p);

            simdovl::compute_prim_pg_overlap_8(buffer, 143, 0, 3, 13, 22, 128, 133, ncols, p);

            simdovl::compute_prim_ph_overlap_5(buffer, 158, 0, 3, 22, 34, 133, 143, ncols, p);

            simdovl::compute_prim_pi_overlap_2(buffer, 181, 0, 3, 34, 51, 143, 158, ncols, p);

            compute_prim_ps_kinetic_energy_0(buffer, 194, 0, 64, 122, ncols, alpha, beta, p);

            compute_prim_pp_kinetic_energy_2(buffer, 197, 3, 64, 125, 194, ncols, alpha, beta, p);

            compute_prim_pd_kinetic_energy_7(buffer, 200, 0, 65, 68, 128, ncols, alpha, beta, p);

            compute_prim_pf_kinetic_energy_10(buffer, 204, 0, 3, 68, 71, 133, 200, ncols, alpha, beta, p);

            compute_prim_pg_kinetic_energy_8(buffer, 209, 0, 3, 71, 80, 143, 204, ncols, alpha, beta, p);

            compute_prim_ph_kinetic_energy_5(buffer, 219, 0, 3, 80, 92, 158, 209, ncols, alpha, beta, p);

            compute_prim_pi_kinetic_energy_2(buffer, 241, 0, 92, 109, 181, ncols, alpha, beta, p);

            simdovl::compute_prim_ds_overlap_2(buffer, 254, 0, 6, 122, ncols, p);

            simdovl::compute_prim_dp_overlap_3(buffer, 257, 0, 3, 122, 125, 254, ncols, p);

            simdovl::compute_prim_dd_overlap_7(buffer, 267, 0, 3, 125, 128, 254, 257, ncols, p);

            simdovl::compute_prim_df_overlap_7(buffer, 284, 0, 3, 128, 133, 257, 267, ncols, p);

            simdovl::compute_prim_dg_overlap_5(buffer, 311, 0, 3, 133, 143, 267, 284, ncols, p);

            simdovl::compute_prim_dh_overlap_3(buffer, 351, 0, 3, 143, 158, 284, 311, ncols, p);

            simdovl::compute_prim_di_overlap_1(buffer, 427, 0, 3, 158, 181, 311, 351, ncols, p);

            compute_prim_ds_kinetic_energy_1(buffer, 503, 0, 6, 64, 194, 254, ncols, alpha, beta, p);

            compute_prim_dp_kinetic_energy_3(buffer, 506, 3, 194, 257, 503, ncols, alpha, beta, p);

            compute_prim_dd_kinetic_energy_7(buffer, 515, 0, 3, 197, 200, 254, 267, 503, 506, ncols, alpha, beta, p);

            compute_prim_df_kinetic_energy_7(buffer, 530, 0, 3, 200, 204, 257, 284, 506, 515, ncols, alpha, beta, p);

            compute_prim_dg_kinetic_energy_5(buffer, 555, 0, 3, 204, 209, 267, 311, 515, 530, ncols, alpha, beta, p);

            compute_prim_dh_kinetic_energy_3(buffer, 595, 0, 3, 209, 219, 284, 351, 530, 555, ncols, alpha, beta, p);

            compute_prim_di_kinetic_energy_1(buffer, 671, 0, 3, 219, 241, 311, 427, 555, 595, ncols, alpha, beta, p);

            simdovl::compute_prim_fs_overlap_6(buffer, 747, 0, 122, 254, ncols, p);

            simdovl::compute_prim_fp_overlap_3(buffer, 753, 0, 3, 254, 257, 747, ncols, p);

            simdovl::compute_prim_fd_overlap_3(buffer, 766, 0, 3, 128, 257, 267, 747, 753, ncols, p);

            simdovl::compute_prim_ff_overlap_3(buffer, 792, 0, 3, 133, 267, 284, 753, 766, ncols, p);

            simdovl::compute_prim_fg_overlap_2(buffer, 840, 0, 3, 143, 284, 311, 766, 792, ncols, p);

            simdovl::compute_prim_fh_overlap_1(buffer, 921, 0, 3, 158, 311, 351, 792, 840, ncols, p);

            simdovl::compute_prim_fi_overlap_0(buffer, 1068, 0, 3, 181, 351, 427, 840, 921, ncols, p);

            compute_prim_fs_kinetic_energy_2(buffer, 1348, 0, 122, 194, 503, 747, ncols, alpha, beta, p);

            compute_prim_fp_kinetic_energy_3(buffer, 1353, 3, 503, 753, 1348, ncols, alpha, beta, p);

            compute_prim_fd_kinetic_energy_3(buffer, 1364, 0, 3, 128, 200, 506, 515, 747, 766, 1348, 1353, ncols, alpha, beta, p);

            compute_prim_ff_kinetic_energy_3(buffer, 1381, 0, 3, 133, 204, 515, 530, 753, 792, 1353, 1364, ncols, alpha, beta, p);

            compute_prim_fg_kinetic_energy_2(buffer, 1411, 0, 3, 143, 209, 530, 555, 766, 840, 1364, 1381, ncols, alpha, beta, p);

            compute_prim_fh_kinetic_energy_1(buffer, 1465, 0, 3, 158, 219, 555, 595, 792, 921, 1381, 1411, ncols, alpha, beta, p);

            compute_prim_fi_kinetic_energy_0(buffer, 1584, 0, 3, 181, 241, 595, 671, 840, 1068, 1411, 1465, ncols, alpha, beta, p);

            simdfunc::contract_primitives(buffer, 1864, 1584, 280, ncols);
        }
    }

    simdtrf::transform_fi(values, nvalues, buffer, 1864, nmax);

    for (size_t m = 0; m < 91; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdkin
