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


#include "SimdNuclearPotentialRecIH.hpp"

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdNuclearPotentialVrrRecDS.hpp"
#include "SimdNuclearPotentialVrrRecFS.hpp"
#include "SimdNuclearPotentialVrrRecGS.hpp"
#include "SimdNuclearPotentialVrrRecHS.hpp"
#include "SimdNuclearPotentialVrrRecIS.hpp"
#include "SimdNuclearPotentialVrrRecKS.hpp"
#include "SimdNuclearPotentialVrrRecLS.hpp"
#include "SimdNuclearPotentialVrrRecMS.hpp"
#include "SimdNuclearPotentialVrrRecNS.hpp"
#include "SimdNuclearPotentialVrrRecOS.hpp"
#include "SimdNuclearPotentialVrrRecPS.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIG.hpp"
#include "SimdTransferIH.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKG.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLF.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMD.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransferNP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_ih_nuclear_potential(double                    *values,
                             const size_t               nvalues,
                             const CBasisFunction      &bra,
                             const CBasisFunction      &ket,
                             const CSimdMatrix         &coordinates,
                             const std::vector<double> &charges,
                             const std::vector<double> &points,
                             CSimdMatrix               &buffer,
                             const double               threshold) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_ih_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_ih_nuclear_potential: Expecting three coordinates for each charge"));

    if (charges.empty())
    {
        std::fill(values, values + 143 * nvalues, 0.0);

        return;
    }

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    // NOTE: the pairs of primitives are screened with the threshold of the
    // integrals divided by their number and by the number of charges, as every
    // integral is a sum over both and the error of a sum is bounded by the
    // number of its terms.

    const auto terms = static_cast<double>(nprims * charges.size());

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_nuclear_potential_primitive_bound, threshold / terms);

    const auto nmax = simdfunc::prepare_buffer(buffer, 6300, 1372, 308, dimensions);

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

            const auto fnpot = 2.0 * mathconst::pi_value() / p * a_norms[i] * b_norms[j];

            const auto fa = -b_exps[j] / p;

            const auto fc = b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto fz = fnpot * charges[ic];

                simdfunc::compute_pc(buffer, coordinates, 3, points, ic, ncols, fc);

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 6, 3, 11, ncols,
                                                          fz, mu, p);

                compute_prim_ps_nuclear_potential_0(buffer, 19, 0, 3, 7, 8, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 22, 0, 3, 8, 9, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 25, 0, 3, 9, 10, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 28, 0, 3, 10, 11, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 31, 0, 3, 11, 12, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 34, 0, 3, 12, 13, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 37, 0, 3, 13, 14, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 40, 0, 3, 14, 15, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 43, 0, 3, 15, 16, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 46, 0, 3, 16, 17, ncols);

                compute_prim_ps_nuclear_potential_0(buffer, 49, 0, 3, 17, 18, ncols);

                compute_prim_ds_nuclear_potential_0(buffer, 52, 0, 3, 7, 8, 19, 22, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 58, 0, 3, 8, 9, 22, 25, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 64, 0, 3, 9, 10, 25, 28, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 70, 0, 3, 10, 11, 28, 31, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 76, 0, 3, 11, 12, 31, 34, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 82, 0, 3, 12, 13, 34, 37, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 88, 0, 3, 13, 14, 37, 40, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 94, 0, 3, 14, 15, 40, 43, ncols, p);

                compute_prim_ds_nuclear_potential_0(buffer, 100, 0, 3, 15, 16, 43, 46, ncols,
                                                    p);

                compute_prim_ds_nuclear_potential_0(buffer, 106, 0, 3, 16, 17, 46, 49, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 112, 0, 3, 19, 22, 52, 58, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 122, 0, 3, 22, 25, 58, 64, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 132, 0, 3, 25, 28, 64, 70, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 142, 0, 3, 28, 31, 70, 76, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 152, 0, 3, 31, 34, 76, 82, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 162, 0, 3, 34, 37, 82, 88, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 172, 0, 3, 37, 40, 88, 94, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 182, 0, 3, 40, 43, 94, 100, ncols,
                                                    p);

                compute_prim_fs_nuclear_potential_0(buffer, 192, 0, 3, 43, 46, 100, 106, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 202, 0, 3, 52, 58, 112, 122, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 217, 0, 3, 58, 64, 122, 132, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 232, 0, 3, 64, 70, 132, 142, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 247, 0, 3, 70, 76, 142, 152, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 262, 0, 3, 76, 82, 152, 162, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 277, 0, 3, 82, 88, 162, 172, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 292, 0, 3, 88, 94, 172, 182, ncols,
                                                    p);

                compute_prim_gs_nuclear_potential_0(buffer, 307, 0, 3, 94, 100, 182, 192, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 322, 0, 3, 112, 122, 202, 217, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 343, 0, 3, 122, 132, 217, 232, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 364, 0, 3, 132, 142, 232, 247, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 385, 0, 3, 142, 152, 247, 262, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 406, 0, 3, 152, 162, 262, 277, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 427, 0, 3, 162, 172, 277, 292, ncols,
                                                    p);

                compute_prim_hs_nuclear_potential_0(buffer, 448, 0, 3, 172, 182, 292, 307, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 469, 0, 3, 202, 217, 322, 343, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 497, 0, 3, 217, 232, 343, 364, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 525, 0, 3, 232, 247, 364, 385, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 553, 0, 3, 247, 262, 385, 406, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 581, 0, 3, 262, 277, 406, 427, ncols,
                                                    p);

                compute_prim_is_nuclear_potential_0(buffer, 609, 0, 3, 277, 292, 427, 448, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 637, 0, 3, 322, 343, 469, 497, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 673, 0, 3, 343, 364, 497, 525, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 709, 0, 3, 364, 385, 525, 553, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 745, 0, 3, 385, 406, 553, 581, ncols,
                                                    p);

                compute_prim_ks_nuclear_potential_0(buffer, 781, 0, 3, 406, 427, 581, 609, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 817, 0, 3, 469, 497, 637, 673, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 862, 0, 3, 497, 525, 673, 709, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 907, 0, 3, 525, 553, 709, 745, ncols,
                                                    p);

                compute_prim_ls_nuclear_potential_0(buffer, 952, 0, 3, 553, 581, 745, 781, ncols,
                                                    p);

                compute_prim_ms_nuclear_potential_0(buffer, 997, 0, 3, 637, 673, 817, 862, ncols,
                                                    p);

                compute_prim_ms_nuclear_potential_0(buffer, 1052, 0, 3, 673, 709, 862, 907,
                                                    ncols, p);

                compute_prim_ms_nuclear_potential_0(buffer, 1107, 0, 3, 709, 745, 907, 952,
                                                    ncols, p);

                compute_prim_ns_nuclear_potential_0(buffer, 1162, 0, 3, 817, 862, 997, 1052,
                                                    ncols, p);

                compute_prim_ns_nuclear_potential_0(buffer, 1228, 0, 3, 862, 907, 1052, 1107,
                                                    ncols, p);

                compute_prim_os_nuclear_potential_0(buffer, 1294, 0, 3, 997, 1052, 1162, 1228,
                                                    ncols, p);

                simdfunc::contract_primitives(buffer, 1372, 469, 28, ncols);

                simdfunc::contract_primitives(buffer, 1400, 637, 36, ncols);

                simdfunc::contract_primitives(buffer, 1436, 817, 45, ncols);

                simdfunc::contract_primitives(buffer, 1481, 997, 55, ncols);

                simdfunc::contract_primitives(buffer, 1536, 1162, 66, ncols);

                simdfunc::contract_primitives(buffer, 1602, 1294, 78, ncols);
            }
        }
    }

    simdtrf::compute_hrr_ip(buffer, coordinates, 1680, 1372, 1400, 1, nmax);

    simdtrf::compute_hrr_kp(buffer, coordinates, 1764, 1400, 1436, 1, nmax);

    simdtrf::compute_hrr_lp(buffer, coordinates, 1872, 1436, 1481, 1, nmax);

    simdtrf::compute_hrr_mp(buffer, coordinates, 2007, 1481, 1536, 1, nmax);

    simdtrf::compute_hrr_np(buffer, coordinates, 2172, 1536, 1602, 1, nmax);

    simdtrf::compute_hrr_id(buffer, coordinates, 2370, 1680, 1764, 1, nmax);

    simdtrf::compute_hrr_kd(buffer, coordinates, 2538, 1764, 1872, 1, nmax);

    simdtrf::compute_hrr_ld(buffer, coordinates, 2754, 1872, 2007, 1, nmax);

    simdtrf::compute_hrr_md(buffer, coordinates, 3024, 2007, 2172, 1, nmax);

    simdtrf::compute_hrr_if(buffer, coordinates, 3354, 2370, 2538, 1, nmax);

    simdtrf::compute_hrr_kf(buffer, coordinates, 3634, 2538, 2754, 1, nmax);

    simdtrf::compute_hrr_lf(buffer, coordinates, 3994, 2754, 3024, 1, nmax);

    simdtrf::compute_hrr_ig(buffer, coordinates, 4444, 3354, 3634, 1, nmax);

    simdtrf::compute_hrr_kg(buffer, coordinates, 4864, 3634, 3994, 1, nmax);

    simdtrf::compute_hrr_ih(buffer, coordinates, 5404, 4444, 4864, 1, nmax);

    simdtrf::transform_h_inner(buffer, 5992, 5404, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 5992, 11, nmax);

    for (size_t m = 0; m < 143; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
