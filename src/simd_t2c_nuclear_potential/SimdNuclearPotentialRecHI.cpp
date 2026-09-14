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


#include "SimdNuclearPotentialRecHI.hpp"

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

#include "SimdNuclearPotentialVrrRecSD.hpp"
#include "SimdNuclearPotentialVrrRecSF.hpp"
#include "SimdNuclearPotentialVrrRecSG.hpp"
#include "SimdNuclearPotentialVrrRecSH.hpp"
#include "SimdNuclearPotentialVrrRecSI.hpp"
#include "SimdNuclearPotentialVrrRecSK.hpp"
#include "SimdNuclearPotentialVrrRecSL.hpp"
#include "SimdNuclearPotentialVrrRecSM.hpp"
#include "SimdNuclearPotentialVrrRecSN.hpp"
#include "SimdNuclearPotentialVrrRecSO.hpp"
#include "SimdNuclearPotentialVrrRecSP.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferDL.hpp"
#include "SimdTransferDM.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferFK.hpp"
#include "SimdTransferFL.hpp"
#include "SimdTransferGI.hpp"
#include "SimdTransferGK.hpp"
#include "SimdTransferHI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransferPM.hpp"
#include "SimdTransferPN.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_hi_nuclear_potential(double                    *values,
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
            false, std::string("compute_hi_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_hi_nuclear_potential: Expecting three coordinates for each charge"));

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

    const auto nmax = simdfunc::prepare_buffer(buffer, 6265, 1372, 308, dimensions);

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

            const auto fb = a_exps[i] / p;

            const auto fc = b_exps[j] / p;

            simdfunc::compute_pb(buffer, coordinates, 0, ncols, fb);

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto fz = fnpot * charges[ic];

                simdfunc::compute_pc(buffer, coordinates, 3, points, ic, ncols, fc);

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 6, 3, 11, ncols,
                                                          fz, mu, p);

                compute_prim_sp_nuclear_potential_0(buffer, 19, 0, 3, 7, 8, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 22, 0, 3, 8, 9, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 25, 0, 3, 9, 10, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 28, 0, 3, 10, 11, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 31, 0, 3, 11, 12, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 34, 0, 3, 12, 13, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 37, 0, 3, 13, 14, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 40, 0, 3, 14, 15, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 43, 0, 3, 15, 16, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 46, 0, 3, 16, 17, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 49, 0, 3, 17, 18, ncols);

                compute_prim_sd_nuclear_potential_0(buffer, 52, 0, 3, 7, 19, 8, 22, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 58, 0, 3, 8, 22, 9, 25, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 64, 0, 3, 9, 25, 10, 28, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 70, 0, 3, 10, 28, 11, 31, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 76, 0, 3, 11, 31, 12, 34, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 82, 0, 3, 12, 34, 13, 37, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 88, 0, 3, 13, 37, 14, 40, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 94, 0, 3, 14, 40, 15, 43, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 100, 0, 3, 15, 43, 16, 46, ncols,
                                                    p);

                compute_prim_sd_nuclear_potential_0(buffer, 106, 0, 3, 16, 46, 17, 49, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 112, 0, 3, 19, 52, 22, 58, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 122, 0, 3, 22, 58, 25, 64, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 132, 0, 3, 25, 64, 28, 70, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 142, 0, 3, 28, 70, 31, 76, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 152, 0, 3, 31, 76, 34, 82, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 162, 0, 3, 34, 82, 37, 88, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 172, 0, 3, 37, 88, 40, 94, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 182, 0, 3, 40, 94, 43, 100, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 192, 0, 3, 43, 100, 46, 106, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 202, 0, 3, 52, 112, 58, 122, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 217, 0, 3, 58, 122, 64, 132, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 232, 0, 3, 64, 132, 70, 142, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 247, 0, 3, 70, 142, 76, 152, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 262, 0, 3, 76, 152, 82, 162, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 277, 0, 3, 82, 162, 88, 172, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 292, 0, 3, 88, 172, 94, 182, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 307, 0, 3, 94, 182, 100, 192, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 322, 0, 3, 112, 202, 122, 217, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 343, 0, 3, 122, 217, 132, 232, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 364, 0, 3, 132, 232, 142, 247, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 385, 0, 3, 142, 247, 152, 262, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 406, 0, 3, 152, 262, 162, 277, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 427, 0, 3, 162, 277, 172, 292, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 448, 0, 3, 172, 292, 182, 307, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 469, 0, 3, 202, 322, 217, 343, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 497, 0, 3, 217, 343, 232, 364, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 525, 0, 3, 232, 364, 247, 385, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 553, 0, 3, 247, 385, 262, 406, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 581, 0, 3, 262, 406, 277, 427, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 609, 0, 3, 277, 427, 292, 448, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 637, 0, 3, 322, 469, 343, 497, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 673, 0, 3, 343, 497, 364, 525, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 709, 0, 3, 364, 525, 385, 553, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 745, 0, 3, 385, 553, 406, 581, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 781, 0, 3, 406, 581, 427, 609, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 817, 0, 3, 469, 637, 497, 673, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 862, 0, 3, 497, 673, 525, 709, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 907, 0, 3, 525, 709, 553, 745, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 952, 0, 3, 553, 745, 581, 781, ncols,
                                                    p);

                compute_prim_sm_nuclear_potential_0(buffer, 997, 0, 3, 637, 817, 673, 862, ncols,
                                                    p);

                compute_prim_sm_nuclear_potential_0(buffer, 1052, 0, 3, 673, 862, 709, 907,
                                                    ncols, p);

                compute_prim_sm_nuclear_potential_0(buffer, 1107, 0, 3, 709, 907, 745, 952,
                                                    ncols, p);

                compute_prim_sn_nuclear_potential_0(buffer, 1162, 0, 3, 817, 997, 862, 1052,
                                                    ncols, p);

                compute_prim_sn_nuclear_potential_0(buffer, 1228, 0, 3, 862, 1052, 907, 1107,
                                                    ncols, p);

                compute_prim_so_nuclear_potential_0(buffer, 1294, 0, 3, 997, 1162, 1052, 1228,
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

    simdtrf::compute_hrr_pi(buffer, coordinates, 1680, 1372, 1400, 1, nmax);

    simdtrf::compute_hrr_pk(buffer, coordinates, 1764, 1400, 1436, 1, nmax);

    simdtrf::compute_hrr_pl(buffer, coordinates, 1872, 1436, 1481, 1, nmax);

    simdtrf::compute_hrr_pm(buffer, coordinates, 2007, 1481, 1536, 1, nmax);

    simdtrf::compute_hrr_pn(buffer, coordinates, 2172, 1536, 1602, 1, nmax);

    simdtrf::compute_hrr_di(buffer, coordinates, 2370, 1680, 1764, 1, nmax);

    simdtrf::compute_hrr_dk(buffer, coordinates, 2538, 1764, 1872, 1, nmax);

    simdtrf::compute_hrr_dl(buffer, coordinates, 2754, 1872, 2007, 1, nmax);

    simdtrf::compute_hrr_dm(buffer, coordinates, 3024, 2007, 2172, 1, nmax);

    simdtrf::compute_hrr_fi(buffer, coordinates, 3354, 2370, 2538, 1, nmax);

    simdtrf::compute_hrr_fk(buffer, coordinates, 3634, 2538, 2754, 1, nmax);

    simdtrf::compute_hrr_fl(buffer, coordinates, 3994, 2754, 3024, 1, nmax);

    simdtrf::compute_hrr_gi(buffer, coordinates, 4444, 3354, 3634, 1, nmax);

    simdtrf::compute_hrr_gk(buffer, coordinates, 4864, 3634, 3994, 1, nmax);

    simdtrf::compute_hrr_hi(buffer, coordinates, 5404, 4444, 4864, 1, nmax);

    simdtrf::transform_i_inner(buffer, 5992, 5404, 21, 1, nmax);

    simdtrf::transform_h_outer(values, nvalues, buffer, 5992, 13, nmax);

    for (size_t m = 0; m < 143; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
