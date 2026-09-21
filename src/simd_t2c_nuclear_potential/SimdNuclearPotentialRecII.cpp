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


#include "SimdNuclearPotentialRecII.hpp"

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
#include "SimdNuclearPotentialVrrRecSQ.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferDL.hpp"
#include "SimdTransferDM.hpp"
#include "SimdTransferDN.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferFK.hpp"
#include "SimdTransferFL.hpp"
#include "SimdTransferFM.hpp"
#include "SimdTransferGI.hpp"
#include "SimdTransferGK.hpp"
#include "SimdTransferGL.hpp"
#include "SimdTransferHI.hpp"
#include "SimdTransferHK.hpp"
#include "SimdTransferII.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransferPM.hpp"
#include "SimdTransferPN.hpp"
#include "SimdTransferPO.hpp"
#include "SimdTransformI.hpp"

namespace simdnpot {  // simdnpot namespace

auto
compute_ii_nuclear_potential(double                    *values,
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
            false, std::string("compute_ii_nuclear_potential: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    errors::assertMsgCritical(
        points.size() == 3 * charges.size(),
        std::string("compute_ii_nuclear_potential: Expecting three coordinates for each charge"));

    if (charges.empty())
    {
        std::fill(values, values + 169 * nvalues, 0.0);

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

    const auto nmax = simdfunc::prepare_buffer(buffer, 10298, 1828, 399, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 169 * nvalues, 0.0);

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

            simdfunc::compute_pair_exponent(buffer, coordinates, 6, ncols, mu);

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto fz = fnpot * charges[ic];

                simdfunc::compute_pc(buffer, coordinates, 3, points, ic, ncols, fc);

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 7, 3, 12, ncols,
                                                          fz, 6, p);

                compute_prim_sp_nuclear_potential_0(buffer, 21, 0, 3, 8, 9, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 24, 0, 3, 9, 10, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 27, 0, 3, 10, 11, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 30, 0, 3, 11, 12, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 33, 0, 3, 12, 13, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 36, 0, 3, 13, 14, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 39, 0, 3, 14, 15, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 42, 0, 3, 15, 16, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 45, 0, 3, 16, 17, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 48, 0, 3, 17, 18, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 51, 0, 3, 18, 19, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 54, 0, 3, 19, 20, ncols);

                compute_prim_sd_nuclear_potential_0(buffer, 57, 0, 3, 8, 21, 9, 24, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 63, 0, 3, 9, 24, 10, 27, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 69, 0, 3, 10, 27, 11, 30, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 75, 0, 3, 11, 30, 12, 33, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 81, 0, 3, 12, 33, 13, 36, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 87, 0, 3, 13, 36, 14, 39, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 93, 0, 3, 14, 39, 15, 42, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 99, 0, 3, 15, 42, 16, 45, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 105, 0, 3, 16, 45, 17, 48, ncols,
                                                    p);

                compute_prim_sd_nuclear_potential_0(buffer, 111, 0, 3, 17, 48, 18, 51, ncols,
                                                    p);

                compute_prim_sd_nuclear_potential_0(buffer, 117, 0, 3, 18, 51, 19, 54, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 123, 0, 3, 21, 57, 24, 63, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 133, 0, 3, 24, 63, 27, 69, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 143, 0, 3, 27, 69, 30, 75, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 153, 0, 3, 30, 75, 33, 81, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 163, 0, 3, 33, 81, 36, 87, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 173, 0, 3, 36, 87, 39, 93, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 183, 0, 3, 39, 93, 42, 99, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 193, 0, 3, 42, 99, 45, 105, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 203, 0, 3, 45, 105, 48, 111, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 213, 0, 3, 48, 111, 51, 117, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 223, 0, 3, 57, 123, 63, 133, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 238, 0, 3, 63, 133, 69, 143, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 253, 0, 3, 69, 143, 75, 153, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 268, 0, 3, 75, 153, 81, 163, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 283, 0, 3, 81, 163, 87, 173, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 298, 0, 3, 87, 173, 93, 183, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 313, 0, 3, 93, 183, 99, 193, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 328, 0, 3, 99, 193, 105, 203, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 343, 0, 3, 105, 203, 111, 213, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 358, 0, 3, 123, 223, 133, 238, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 379, 0, 3, 133, 238, 143, 253, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 400, 0, 3, 143, 253, 153, 268, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 421, 0, 3, 153, 268, 163, 283, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 442, 0, 3, 163, 283, 173, 298, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 463, 0, 3, 173, 298, 183, 313, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 484, 0, 3, 183, 313, 193, 328, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 505, 0, 3, 193, 328, 203, 343, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 526, 0, 3, 223, 358, 238, 379, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 554, 0, 3, 238, 379, 253, 400, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 582, 0, 3, 253, 400, 268, 421, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 610, 0, 3, 268, 421, 283, 442, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 638, 0, 3, 283, 442, 298, 463, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 666, 0, 3, 298, 463, 313, 484, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 694, 0, 3, 313, 484, 328, 505, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 722, 0, 3, 358, 526, 379, 554, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 758, 0, 3, 379, 554, 400, 582, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 794, 0, 3, 400, 582, 421, 610, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 830, 0, 3, 421, 610, 442, 638, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 866, 0, 3, 442, 638, 463, 666, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 902, 0, 3, 463, 666, 484, 694, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 938, 0, 3, 526, 722, 554, 758, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 983, 0, 3, 554, 758, 582, 794, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 1028, 0, 3, 582, 794, 610, 830,
                                                    ncols, p);

                compute_prim_sl_nuclear_potential_0(buffer, 1073, 0, 3, 610, 830, 638, 866,
                                                    ncols, p);

                compute_prim_sl_nuclear_potential_0(buffer, 1118, 0, 3, 638, 866, 666, 902,
                                                    ncols, p);

                compute_prim_sm_nuclear_potential_0(buffer, 1163, 0, 3, 722, 938, 758, 983,
                                                    ncols, p);

                compute_prim_sm_nuclear_potential_0(buffer, 1218, 0, 3, 758, 983, 794, 1028,
                                                    ncols, p);

                compute_prim_sm_nuclear_potential_0(buffer, 1273, 0, 3, 794, 1028, 830, 1073,
                                                    ncols, p);

                compute_prim_sm_nuclear_potential_0(buffer, 1328, 0, 3, 830, 1073, 866, 1118,
                                                    ncols, p);

                compute_prim_sn_nuclear_potential_0(buffer, 1383, 0, 3, 938, 1163, 983, 1218,
                                                    ncols, p);

                compute_prim_sn_nuclear_potential_0(buffer, 1449, 0, 3, 983, 1218, 1028, 1273,
                                                    ncols, p);

                compute_prim_sn_nuclear_potential_0(buffer, 1515, 0, 3, 1028, 1273, 1073, 1328,
                                                    ncols, p);

                compute_prim_so_nuclear_potential_0(buffer, 1581, 0, 3, 1163, 1383, 1218, 1449,
                                                    ncols, p);

                compute_prim_so_nuclear_potential_0(buffer, 1659, 0, 3, 1218, 1449, 1273, 1515,
                                                    ncols, p);

                compute_prim_sq_nuclear_potential_0(buffer, 1737, 0, 3, 1383, 1581, 1449, 1659,
                                                    ncols, p);

                simdfunc::contract_primitives(buffer, 1828, 526, 28, ncols);

                simdfunc::contract_primitives(buffer, 1856, 722, 36, ncols);

                simdfunc::contract_primitives(buffer, 1892, 938, 45, ncols);

                simdfunc::contract_primitives(buffer, 1937, 1163, 55, ncols);

                simdfunc::contract_primitives(buffer, 1992, 1383, 66, ncols);

                simdfunc::contract_primitives(buffer, 2058, 1581, 78, ncols);

                simdfunc::contract_primitives(buffer, 2136, 1737, 91, ncols);
            }
        }
    }

    simdtrf::compute_hrr_pi(buffer, coordinates, 2227, 1828, 1856, 1, nmax);

    simdtrf::compute_hrr_pk(buffer, coordinates, 2311, 1856, 1892, 1, nmax);

    simdtrf::compute_hrr_pl(buffer, coordinates, 2419, 1892, 1937, 1, nmax);

    simdtrf::compute_hrr_pm(buffer, coordinates, 2554, 1937, 1992, 1, nmax);

    simdtrf::compute_hrr_pn(buffer, coordinates, 2719, 1992, 2058, 1, nmax);

    simdtrf::compute_hrr_po(buffer, coordinates, 2917, 2058, 2136, 1, nmax);

    simdtrf::compute_hrr_di(buffer, coordinates, 3151, 2227, 2311, 1, nmax);

    simdtrf::compute_hrr_dk(buffer, coordinates, 3319, 2311, 2419, 1, nmax);

    simdtrf::compute_hrr_dl(buffer, coordinates, 3535, 2419, 2554, 1, nmax);

    simdtrf::compute_hrr_dm(buffer, coordinates, 3805, 2554, 2719, 1, nmax);

    simdtrf::compute_hrr_dn(buffer, coordinates, 4135, 2719, 2917, 1, nmax);

    simdtrf::compute_hrr_fi(buffer, coordinates, 4531, 3151, 3319, 1, nmax);

    simdtrf::compute_hrr_fk(buffer, coordinates, 4811, 3319, 3535, 1, nmax);

    simdtrf::compute_hrr_fl(buffer, coordinates, 5171, 3535, 3805, 1, nmax);

    simdtrf::compute_hrr_fm(buffer, coordinates, 5621, 3805, 4135, 1, nmax);

    simdtrf::compute_hrr_gi(buffer, coordinates, 6171, 4531, 4811, 1, nmax);

    simdtrf::compute_hrr_gk(buffer, coordinates, 6591, 4811, 5171, 1, nmax);

    simdtrf::compute_hrr_gl(buffer, coordinates, 7131, 5171, 5621, 1, nmax);

    simdtrf::compute_hrr_hi(buffer, coordinates, 7806, 6171, 6591, 1, nmax);

    simdtrf::compute_hrr_hk(buffer, coordinates, 8394, 6591, 7131, 1, nmax);

    simdtrf::compute_hrr_ii(buffer, coordinates, 9150, 7806, 8394, 1, nmax);

    simdtrf::transform_i_inner(buffer, 9934, 9150, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 9934, 13, nmax);

    for (size_t m = 0; m < 169; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
