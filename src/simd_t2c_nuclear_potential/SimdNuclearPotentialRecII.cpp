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

    const auto nmax = simdfunc::prepare_buffer(buffer, 10297, 1827, 399, dimensions);

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

            for (size_t ic = 0; ic < charges.size(); ic++)
            {
                const auto fz = fnpot * charges[ic];

                simdfunc::compute_pc(buffer, coordinates, 3, points, ic, ncols, fc);

                simdfunc::compute_full_npot_boys_function(buffer, coordinates, 6, 3, 12, ncols,
                                                          fz, mu, p);

                compute_prim_sp_nuclear_potential_0(buffer, 20, 0, 3, 7, 8, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 23, 0, 3, 8, 9, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 26, 0, 3, 9, 10, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 29, 0, 3, 10, 11, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 32, 0, 3, 11, 12, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 35, 0, 3, 12, 13, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 38, 0, 3, 13, 14, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 41, 0, 3, 14, 15, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 44, 0, 3, 15, 16, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 47, 0, 3, 16, 17, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 50, 0, 3, 17, 18, ncols);

                compute_prim_sp_nuclear_potential_0(buffer, 53, 0, 3, 18, 19, ncols);

                compute_prim_sd_nuclear_potential_0(buffer, 56, 0, 3, 7, 20, 8, 23, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 62, 0, 3, 8, 23, 9, 26, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 68, 0, 3, 9, 26, 10, 29, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 74, 0, 3, 10, 29, 11, 32, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 80, 0, 3, 11, 32, 12, 35, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 86, 0, 3, 12, 35, 13, 38, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 92, 0, 3, 13, 38, 14, 41, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 98, 0, 3, 14, 41, 15, 44, ncols, p);

                compute_prim_sd_nuclear_potential_0(buffer, 104, 0, 3, 15, 44, 16, 47, ncols,
                                                    p);

                compute_prim_sd_nuclear_potential_0(buffer, 110, 0, 3, 16, 47, 17, 50, ncols,
                                                    p);

                compute_prim_sd_nuclear_potential_0(buffer, 116, 0, 3, 17, 50, 18, 53, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 122, 0, 3, 20, 56, 23, 62, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 132, 0, 3, 23, 62, 26, 68, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 142, 0, 3, 26, 68, 29, 74, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 152, 0, 3, 29, 74, 32, 80, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 162, 0, 3, 32, 80, 35, 86, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 172, 0, 3, 35, 86, 38, 92, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 182, 0, 3, 38, 92, 41, 98, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 192, 0, 3, 41, 98, 44, 104, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 202, 0, 3, 44, 104, 47, 110, ncols,
                                                    p);

                compute_prim_sf_nuclear_potential_0(buffer, 212, 0, 3, 47, 110, 50, 116, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 222, 0, 3, 56, 122, 62, 132, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 237, 0, 3, 62, 132, 68, 142, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 252, 0, 3, 68, 142, 74, 152, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 267, 0, 3, 74, 152, 80, 162, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 282, 0, 3, 80, 162, 86, 172, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 297, 0, 3, 86, 172, 92, 182, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 312, 0, 3, 92, 182, 98, 192, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 327, 0, 3, 98, 192, 104, 202, ncols,
                                                    p);

                compute_prim_sg_nuclear_potential_0(buffer, 342, 0, 3, 104, 202, 110, 212, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 357, 0, 3, 122, 222, 132, 237, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 378, 0, 3, 132, 237, 142, 252, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 399, 0, 3, 142, 252, 152, 267, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 420, 0, 3, 152, 267, 162, 282, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 441, 0, 3, 162, 282, 172, 297, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 462, 0, 3, 172, 297, 182, 312, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 483, 0, 3, 182, 312, 192, 327, ncols,
                                                    p);

                compute_prim_sh_nuclear_potential_0(buffer, 504, 0, 3, 192, 327, 202, 342, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 525, 0, 3, 222, 357, 237, 378, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 553, 0, 3, 237, 378, 252, 399, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 581, 0, 3, 252, 399, 267, 420, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 609, 0, 3, 267, 420, 282, 441, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 637, 0, 3, 282, 441, 297, 462, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 665, 0, 3, 297, 462, 312, 483, ncols,
                                                    p);

                compute_prim_si_nuclear_potential_0(buffer, 693, 0, 3, 312, 483, 327, 504, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 721, 0, 3, 357, 525, 378, 553, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 757, 0, 3, 378, 553, 399, 581, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 793, 0, 3, 399, 581, 420, 609, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 829, 0, 3, 420, 609, 441, 637, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 865, 0, 3, 441, 637, 462, 665, ncols,
                                                    p);

                compute_prim_sk_nuclear_potential_0(buffer, 901, 0, 3, 462, 665, 483, 693, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 937, 0, 3, 525, 721, 553, 757, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 982, 0, 3, 553, 757, 581, 793, ncols,
                                                    p);

                compute_prim_sl_nuclear_potential_0(buffer, 1027, 0, 3, 581, 793, 609, 829,
                                                    ncols, p);

                compute_prim_sl_nuclear_potential_0(buffer, 1072, 0, 3, 609, 829, 637, 865,
                                                    ncols, p);

                compute_prim_sl_nuclear_potential_0(buffer, 1117, 0, 3, 637, 865, 665, 901,
                                                    ncols, p);

                compute_prim_sm_nuclear_potential_0(buffer, 1162, 0, 3, 721, 937, 757, 982,
                                                    ncols, p);

                compute_prim_sm_nuclear_potential_0(buffer, 1217, 0, 3, 757, 982, 793, 1027,
                                                    ncols, p);

                compute_prim_sm_nuclear_potential_0(buffer, 1272, 0, 3, 793, 1027, 829, 1072,
                                                    ncols, p);

                compute_prim_sm_nuclear_potential_0(buffer, 1327, 0, 3, 829, 1072, 865, 1117,
                                                    ncols, p);

                compute_prim_sn_nuclear_potential_0(buffer, 1382, 0, 3, 937, 1162, 982, 1217,
                                                    ncols, p);

                compute_prim_sn_nuclear_potential_0(buffer, 1448, 0, 3, 982, 1217, 1027, 1272,
                                                    ncols, p);

                compute_prim_sn_nuclear_potential_0(buffer, 1514, 0, 3, 1027, 1272, 1072, 1327,
                                                    ncols, p);

                compute_prim_so_nuclear_potential_0(buffer, 1580, 0, 3, 1162, 1382, 1217, 1448,
                                                    ncols, p);

                compute_prim_so_nuclear_potential_0(buffer, 1658, 0, 3, 1217, 1448, 1272, 1514,
                                                    ncols, p);

                compute_prim_sq_nuclear_potential_0(buffer, 1736, 0, 3, 1382, 1580, 1448, 1658,
                                                    ncols, p);

                simdfunc::contract_primitives(buffer, 1827, 525, 28, ncols);

                simdfunc::contract_primitives(buffer, 1855, 721, 36, ncols);

                simdfunc::contract_primitives(buffer, 1891, 937, 45, ncols);

                simdfunc::contract_primitives(buffer, 1936, 1162, 55, ncols);

                simdfunc::contract_primitives(buffer, 1991, 1382, 66, ncols);

                simdfunc::contract_primitives(buffer, 2057, 1580, 78, ncols);

                simdfunc::contract_primitives(buffer, 2135, 1736, 91, ncols);
            }
        }
    }

    simdtrf::compute_hrr_pi(buffer, coordinates, 2226, 1827, 1855, 1, nmax);

    simdtrf::compute_hrr_pk(buffer, coordinates, 2310, 1855, 1891, 1, nmax);

    simdtrf::compute_hrr_pl(buffer, coordinates, 2418, 1891, 1936, 1, nmax);

    simdtrf::compute_hrr_pm(buffer, coordinates, 2553, 1936, 1991, 1, nmax);

    simdtrf::compute_hrr_pn(buffer, coordinates, 2718, 1991, 2057, 1, nmax);

    simdtrf::compute_hrr_po(buffer, coordinates, 2916, 2057, 2135, 1, nmax);

    simdtrf::compute_hrr_di(buffer, coordinates, 3150, 2226, 2310, 1, nmax);

    simdtrf::compute_hrr_dk(buffer, coordinates, 3318, 2310, 2418, 1, nmax);

    simdtrf::compute_hrr_dl(buffer, coordinates, 3534, 2418, 2553, 1, nmax);

    simdtrf::compute_hrr_dm(buffer, coordinates, 3804, 2553, 2718, 1, nmax);

    simdtrf::compute_hrr_dn(buffer, coordinates, 4134, 2718, 2916, 1, nmax);

    simdtrf::compute_hrr_fi(buffer, coordinates, 4530, 3150, 3318, 1, nmax);

    simdtrf::compute_hrr_fk(buffer, coordinates, 4810, 3318, 3534, 1, nmax);

    simdtrf::compute_hrr_fl(buffer, coordinates, 5170, 3534, 3804, 1, nmax);

    simdtrf::compute_hrr_fm(buffer, coordinates, 5620, 3804, 4134, 1, nmax);

    simdtrf::compute_hrr_gi(buffer, coordinates, 6170, 4530, 4810, 1, nmax);

    simdtrf::compute_hrr_gk(buffer, coordinates, 6590, 4810, 5170, 1, nmax);

    simdtrf::compute_hrr_gl(buffer, coordinates, 7130, 5170, 5620, 1, nmax);

    simdtrf::compute_hrr_hi(buffer, coordinates, 7805, 6170, 6590, 1, nmax);

    simdtrf::compute_hrr_hk(buffer, coordinates, 8393, 6590, 7130, 1, nmax);

    simdtrf::compute_hrr_ii(buffer, coordinates, 9149, 7805, 8393, 1, nmax);

    simdtrf::transform_i_inner(buffer, 9933, 9149, 28, 1, nmax);

    simdtrf::transform_i_outer(values, nvalues, buffer, 9933, 13, nmax);

    for (size_t m = 0; m < 169; m++)
    {
        auto *pv = values + m * nvalues;

        std::fill(pv + nmax, pv + nvalues, 0.0);
    }
}

}  // namespace simdnpot
