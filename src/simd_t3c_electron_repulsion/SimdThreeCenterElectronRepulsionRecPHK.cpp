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


#include "SimdThreeCenterElectronRepulsionRecPHK.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSII.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_phk_three_center_electron_repulsion(double               *values,
                                            const size_t          npairs,
                                            const size_t          natoms,
                                            const CBasisFunction &a_function,
                                            const CBasisFunction &b_function,
                                            const CBasisFunction &c_function,
                                            const CSimdMatrix    &coordinates,
                                            const CSimdMatrix    &c_coordinates,
                                            CSimdMatrix          &buffer,
                                            const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_phk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (npairs == 0 || natoms == 0) return;

    const auto &a_exps = a_function.exponents();

    const auto &b_exps = b_function.exponents();

    const auto &c_exps = c_function.exponents();

    const auto &a_norms = a_function.normalization_factors();

    const auto &b_norms = b_function.normalization_factors();

    const auto &c_norms = c_function.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprim_c = c_exps.size();

    const auto nprims = nprim_a * nprim_b * nprim_c;

    // NOTE: the bound neglects the position of the atom on the ket side, so the
    // columns that survive are the same for every one of them and are counted
    // once here rather than inside the loop over them.

    const auto dimensions = simdfunc::make_column_dimensions(
        a_function, b_function, c_function, npairs, coordinates,
        screenfunc::three_center_electron_repulsion_primitive_bound,
        threshold / static_cast<double>(nprims));

    const auto nmax = simdfunc::prepare_buffer(buffer, 42569, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 495 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 42569, 38630, 2079, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pb(buffer, coordinates, 0, nmax, fb);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 3, n, nmax, fc);

                simdfunc::compute_pair_exponent(buffer, coordinates, 6, nmax, mu);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 7, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13}, ncols,
                                                        fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 57, 0, 3, 8, 9,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 63, 0, 3, 9, 10,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 69, 0, 3, 10, 11,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 75, 0, 3, 11, 12,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 81, 0, 3, 12, 13,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 87, 0, 3, 13, 14,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 93, 0, 3, 14, 15,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 99, 0, 3, 15, 16,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 105, 0, 3, 16, 17,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 111, 0, 3, 17, 18,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 117, 0, 3, 18, 19,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 123, 0, 3, 21, 24,
                                                                       57, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 133, 0, 3, 24, 27,
                                                                       63, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 143, 0, 3, 27, 30,
                                                                       69, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 153, 0, 3, 30, 33,
                                                                       75, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 163, 0, 3, 33, 36,
                                                                       81, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 173, 0, 3, 36, 39,
                                                                       87, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 183, 0, 3, 39, 42,
                                                                       93, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 193, 0, 3, 42, 45,
                                                                       99, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 45, 48,
                                                                       105, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 48, 51,
                                                                       111, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 223, 0, 3, 57, 63,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 63, 69,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 253, 0, 3, 69, 75,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 75, 81,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 283, 0, 3, 81, 87,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 87, 93,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 313, 0, 3, 93, 99,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 99,
                                                                       105, 193, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 343, 0, 3, 105,
                                                                       111, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 123,
                                                                       133, 223, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 379, 0, 3, 133,
                                                                       143, 238, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 400, 0, 3, 143,
                                                                       153, 253, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 421, 0, 3, 153,
                                                                       163, 268, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 442, 0, 3, 163,
                                                                       173, 283, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 463, 0, 3, 173,
                                                                       183, 298, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 484, 0, 3, 183,
                                                                       193, 313, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 505, 0, 3, 193,
                                                                       203, 328, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 526, 0, 3, 223,
                                                                       238, 358, 379, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 554, 0, 3, 238,
                                                                       253, 379, 400, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 582, 0, 3, 253,
                                                                       268, 400, 421, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 610, 0, 3, 268,
                                                                       283, 421, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 638, 0, 3, 283,
                                                                       298, 442, 463, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 666, 0, 3, 298,
                                                                       313, 463, 484, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 694, 0, 3, 313,
                                                                       328, 484, 505, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 722, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 725, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 728, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 731, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 734, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 737, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 740, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 743, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 746, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 749, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 752, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 755, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 758, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 761, 3, 10, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 770, 3, 11, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 779, 3, 12, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 788, 3, 13, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 797, 3, 14, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 806, 3, 15, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 815, 3, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 824, 3, 17, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 833, 3, 18, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 842, 3, 19, 54,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 851, 3, 21, 57,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 869, 3, 24, 63,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 887, 3, 27, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 905, 3, 30, 75,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 923, 3, 33, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 941, 3, 36, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 959, 3, 39, 93,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 977, 3, 42, 99,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 995, 3, 45, 105,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1013, 3, 48, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1031, 3, 51, 117,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1049, 3, 57, 123,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1079, 3, 63, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1109, 3, 69, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1139, 3, 75, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1169, 3, 81, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1199, 3, 87, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1229, 3, 93, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1259, 3, 99, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1289, 3, 105, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1319, 3, 111, 213,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1349, 3, 123, 223,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1394, 3, 133, 238,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1439, 3, 143, 253,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1484, 3, 153, 268,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1529, 3, 163, 283,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1574, 3, 173, 298,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1619, 3, 183, 313,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1664, 3, 193, 328,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1709, 3, 203, 343,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1754, 3, 223, 358,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1817, 3, 238, 379,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1880, 3, 253, 400,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1943, 3, 268, 421,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2006, 3, 283, 442,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2069, 3, 298, 463,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2132, 3, 313, 484,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2195, 3, 328, 505,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2258, 3, 358, 526,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2342, 3, 379, 554,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2426, 3, 400, 582,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2510, 3, 421, 610,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2594, 3, 442, 638,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2678, 3, 463, 666,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2762, 3, 484, 694,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2846, 3, 8, 9,
                                                                       728, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2852, 3, 9, 10,
                                                                       731, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2858, 3, 10, 11,
                                                                       734, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2864, 3, 11, 12,
                                                                       737, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2870, 3, 12, 13,
                                                                       740, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2876, 3, 13, 14,
                                                                       743, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2882, 3, 14, 15,
                                                                       746, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2888, 3, 15, 16,
                                                                       749, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2894, 3, 16, 17,
                                                                       752, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2900, 3, 17, 18,
                                                                       755, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2906, 3, 18, 19,
                                                                       758, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2912, 0, 3, 2846,
                                                                       728, 2852, 761, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2930, 0, 3, 2852,
                                                                       731, 2858, 770, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2948, 0, 3, 2858,
                                                                       734, 2864, 779, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2966, 0, 3, 2864,
                                                                       737, 2870, 788, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2984, 0, 3, 2870,
                                                                       740, 2876, 797, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3002, 0, 3, 2876,
                                                                       743, 2882, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3020, 0, 3, 2882,
                                                                       746, 2888, 815, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3038, 0, 3, 2888,
                                                                       749, 2894, 824, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3056, 0, 3, 2894,
                                                                       752, 2900, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3074, 0, 3, 2900,
                                                                       755, 2906, 842, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3092, 0, 3, 2912,
                                                                       761, 2930, 57, 63, 887,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3128, 0, 3, 2930,
                                                                       770, 2948, 63, 69, 905,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3164, 0, 3, 2948,
                                                                       779, 2966, 69, 75, 923,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3200, 0, 3, 2966,
                                                                       788, 2984, 75, 81, 941,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3236, 0, 3, 2984,
                                                                       797, 3002, 81, 87, 959,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3272, 0, 3, 3002,
                                                                       806, 3020, 87, 93, 977,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3308, 0, 3, 3020,
                                                                       815, 3038, 93, 99, 995,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3344, 0, 3, 3038,
                                                                       824, 3056, 99, 105, 1013,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3380, 0, 3, 3056,
                                                                       833, 3074, 105, 111, 1031,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3416, 0, 3, 3092,
                                                                       887, 3128, 123, 133, 1109,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3476, 0, 3, 3128,
                                                                       905, 3164, 133, 143, 1139,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3536, 0, 3, 3164,
                                                                       923, 3200, 143, 153, 1169,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3596, 0, 3, 3200,
                                                                       941, 3236, 153, 163, 1199,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3656, 0, 3, 3236,
                                                                       959, 3272, 163, 173, 1229,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3716, 0, 3, 3272,
                                                                       977, 3308, 173, 183, 1259,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3776, 0, 3, 3308,
                                                                       995, 3344, 183, 193, 1289,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3836, 0, 3, 3344,
                                                                       1013, 3380, 193, 203,
                                                                       1319, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3896, 0, 3, 3416,
                                                                       1109, 3476, 223, 238,
                                                                       1439, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3986, 0, 3, 3476,
                                                                       1139, 3536, 238, 253,
                                                                       1484, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4076, 0, 3, 3536,
                                                                       1169, 3596, 253, 268,
                                                                       1529, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4166, 0, 3, 3596,
                                                                       1199, 3656, 268, 283,
                                                                       1574, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4256, 0, 3, 3656,
                                                                       1229, 3716, 283, 298,
                                                                       1619, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4346, 0, 3, 3716,
                                                                       1259, 3776, 298, 313,
                                                                       1664, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4436, 0, 3, 3776,
                                                                       1289, 3836, 313, 328,
                                                                       1709, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4526, 0, 3, 3896,
                                                                       1439, 3986, 358, 379,
                                                                       1880, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4652, 0, 3, 3986,
                                                                       1484, 4076, 379, 400,
                                                                       1943, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4778, 0, 3, 4076,
                                                                       1529, 4166, 400, 421,
                                                                       2006, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 4904, 0, 3, 4166,
                                                                       1574, 4256, 421, 442,
                                                                       2069, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5030, 0, 3, 4256,
                                                                       1619, 4346, 442, 463,
                                                                       2132, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5156, 0, 3, 4346,
                                                                       1664, 4436, 463, 484,
                                                                       2195, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 5282, 0, 3, 4526,
                                                                       1880, 4652, 526, 554,
                                                                       2426, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 5450, 0, 3, 4652,
                                                                       1943, 4778, 554, 582,
                                                                       2510, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 5618, 0, 3, 4778,
                                                                       2006, 4904, 582, 610,
                                                                       2594, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 5786, 0, 3, 4904,
                                                                       2069, 5030, 610, 638,
                                                                       2678, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 5954, 0, 3, 5030,
                                                                       2132, 5156, 638, 666,
                                                                       2762, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6122, 3, 722, 725,
                                                                       2846, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6132, 3, 725, 728,
                                                                       2852, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6142, 3, 728, 731,
                                                                       2858, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6152, 3, 731, 734,
                                                                       2864, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6162, 3, 734, 737,
                                                                       2870, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6172, 3, 737, 740,
                                                                       2876, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6182, 3, 740, 743,
                                                                       2882, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6192, 3, 743, 746,
                                                                       2888, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6202, 3, 746, 749,
                                                                       2894, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6212, 3, 749, 752,
                                                                       2900, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6222, 3, 752, 755,
                                                                       2906, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6232, 0, 3, 6122,
                                                                       2846, 6132, 2912, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6262, 0, 3, 6132,
                                                                       2852, 6142, 2930, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6292, 0, 3, 6142,
                                                                       2858, 6152, 2948, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6322, 0, 3, 6152,
                                                                       2864, 6162, 2966, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6352, 0, 3, 6162,
                                                                       2870, 6172, 2984, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6382, 0, 3, 6172,
                                                                       2876, 6182, 3002, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6412, 0, 3, 6182,
                                                                       2882, 6192, 3020, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6442, 0, 3, 6192,
                                                                       2888, 6202, 3038, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6472, 0, 3, 6202,
                                                                       2894, 6212, 3056, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6502, 0, 3, 6212,
                                                                       2900, 6222, 3074, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6532, 0, 3, 6232,
                                                                       2912, 6262, 851, 869,
                                                                       3092, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6592, 0, 3, 6262,
                                                                       2930, 6292, 869, 887,
                                                                       3128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6652, 0, 3, 6292,
                                                                       2948, 6322, 887, 905,
                                                                       3164, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6712, 0, 3, 6322,
                                                                       2966, 6352, 905, 923,
                                                                       3200, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6772, 0, 3, 6352,
                                                                       2984, 6382, 923, 941,
                                                                       3236, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6832, 0, 3, 6382,
                                                                       3002, 6412, 941, 959,
                                                                       3272, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6892, 0, 3, 6412,
                                                                       3020, 6442, 959, 977,
                                                                       3308, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 6952, 0, 3, 6442,
                                                                       3038, 6472, 977, 995,
                                                                       3344, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7012, 0, 3, 6472,
                                                                       3056, 6502, 995, 1013,
                                                                       3380, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7072, 0, 3, 6532,
                                                                       3092, 6592, 1049, 1079,
                                                                       3416, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7172, 0, 3, 6592,
                                                                       3128, 6652, 1079, 1109,
                                                                       3476, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7272, 0, 3, 6652,
                                                                       3164, 6712, 1109, 1139,
                                                                       3536, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7372, 0, 3, 6712,
                                                                       3200, 6772, 1139, 1169,
                                                                       3596, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7472, 0, 3, 6772,
                                                                       3236, 6832, 1169, 1199,
                                                                       3656, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7572, 0, 3, 6832,
                                                                       3272, 6892, 1199, 1229,
                                                                       3716, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7672, 0, 3, 6892,
                                                                       3308, 6952, 1229, 1259,
                                                                       3776, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7772, 0, 3, 6952,
                                                                       3344, 7012, 1259, 1289,
                                                                       3836, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 7872, 0, 3, 7072,
                                                                       3416, 7172, 1349, 1394,
                                                                       3896, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8022, 0, 3, 7172,
                                                                       3476, 7272, 1394, 1439,
                                                                       3986, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8172, 0, 3, 7272,
                                                                       3536, 7372, 1439, 1484,
                                                                       4076, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8322, 0, 3, 7372,
                                                                       3596, 7472, 1484, 1529,
                                                                       4166, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8472, 0, 3, 7472,
                                                                       3656, 7572, 1529, 1574,
                                                                       4256, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8622, 0, 3, 7572,
                                                                       3716, 7672, 1574, 1619,
                                                                       4346, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8772, 0, 3, 7672,
                                                                       3776, 7772, 1619, 1664,
                                                                       4436, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 8922, 0, 3, 7872,
                                                                       3896, 8022, 1754, 1817,
                                                                       4526, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 9132, 0, 3, 8022,
                                                                       3986, 8172, 1817, 1880,
                                                                       4652, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 9342, 0, 3, 8172,
                                                                       4076, 8322, 1880, 1943,
                                                                       4778, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 9552, 0, 3, 8322,
                                                                       4166, 8472, 1943, 2006,
                                                                       4904, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 9762, 0, 3, 8472,
                                                                       4256, 8622, 2006, 2069,
                                                                       5030, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 9972, 0, 3, 8622,
                                                                       4346, 8772, 2069, 2132,
                                                                       5156, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 10182, 0, 3, 8922,
                                                                       4526, 9132, 2258, 2342,
                                                                       5282, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 10462, 0, 3, 9132,
                                                                       4652, 9342, 2342, 2426,
                                                                       5450, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 10742, 0, 3, 9342,
                                                                       4778, 9552, 2426, 2510,
                                                                       5618, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 11022, 0, 3, 9552,
                                                                       4904, 9762, 2510, 2594,
                                                                       5786, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 11302, 0, 3, 9762,
                                                                       5030, 9972, 2594, 2678,
                                                                       5954, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11582, 3, 2846,
                                                                       2852, 6142, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11597, 3, 2852,
                                                                       2858, 6152, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11612, 3, 2858,
                                                                       2864, 6162, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11627, 3, 2864,
                                                                       2870, 6172, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11642, 3, 2870,
                                                                       2876, 6182, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11657, 3, 2876,
                                                                       2882, 6192, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11672, 3, 2882,
                                                                       2888, 6202, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11687, 3, 2888,
                                                                       2894, 6212, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 11702, 3, 2894,
                                                                       2900, 6222, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 11717, 0, 3,
                                                                       11582, 6142, 11597, 2912,
                                                                       2930, 6292, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 11762, 0, 3,
                                                                       11597, 6152, 11612, 2930,
                                                                       2948, 6322, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 11807, 0, 3,
                                                                       11612, 6162, 11627, 2948,
                                                                       2966, 6352, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 11852, 0, 3,
                                                                       11627, 6172, 11642, 2966,
                                                                       2984, 6382, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 11897, 0, 3,
                                                                       11642, 6182, 11657, 2984,
                                                                       3002, 6412, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 11942, 0, 3,
                                                                       11657, 6192, 11672, 3002,
                                                                       3020, 6442, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 11987, 0, 3,
                                                                       11672, 6202, 11687, 3020,
                                                                       3038, 6472, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 12032, 0, 3,
                                                                       11687, 6212, 11702, 3038,
                                                                       3056, 6502, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 12077, 0, 3,
                                                                       11717, 6292, 11762, 3092,
                                                                       3128, 6652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 12167, 0, 3,
                                                                       11762, 6322, 11807, 3128,
                                                                       3164, 6712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 12257, 0, 3,
                                                                       11807, 6352, 11852, 3164,
                                                                       3200, 6772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 12347, 0, 3,
                                                                       11852, 6382, 11897, 3200,
                                                                       3236, 6832, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 12437, 0, 3,
                                                                       11897, 6412, 11942, 3236,
                                                                       3272, 6892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 12527, 0, 3,
                                                                       11942, 6442, 11987, 3272,
                                                                       3308, 6952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 12617, 0, 3,
                                                                       11987, 6472, 12032, 3308,
                                                                       3344, 7012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 12707, 0, 3,
                                                                       12077, 6652, 12167, 3416,
                                                                       3476, 7272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 12857, 0, 3,
                                                                       12167, 6712, 12257, 3476,
                                                                       3536, 7372, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 13007, 0, 3,
                                                                       12257, 6772, 12347, 3536,
                                                                       3596, 7472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 13157, 0, 3,
                                                                       12347, 6832, 12437, 3596,
                                                                       3656, 7572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 13307, 0, 3,
                                                                       12437, 6892, 12527, 3656,
                                                                       3716, 7672, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 13457, 0, 3,
                                                                       12527, 6952, 12617, 3716,
                                                                       3776, 7772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 13607, 0, 3,
                                                                       12707, 7272, 12857, 3896,
                                                                       3986, 8172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 13832, 0, 3,
                                                                       12857, 7372, 13007, 3986,
                                                                       4076, 8322, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 14057, 0, 3,
                                                                       13007, 7472, 13157, 4076,
                                                                       4166, 8472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 14282, 0, 3,
                                                                       13157, 7572, 13307, 4166,
                                                                       4256, 8622, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 14507, 0, 3,
                                                                       13307, 7672, 13457, 4256,
                                                                       4346, 8772, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 14732, 0, 3,
                                                                       13607, 8172, 13832, 4526,
                                                                       4652, 9342, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 15047, 0, 3,
                                                                       13832, 8322, 14057, 4652,
                                                                       4778, 9552, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 15362, 0, 3,
                                                                       14057, 8472, 14282, 4778,
                                                                       4904, 9762, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 15677, 0, 3,
                                                                       14282, 8622, 14507, 4904,
                                                                       5030, 9972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 15992, 0, 3,
                                                                       14732, 9342, 15047, 5282,
                                                                       5450, 10742, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 16412, 0, 3,
                                                                       15047, 9552, 15362, 5450,
                                                                       5618, 11022, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 16832, 0, 3,
                                                                       15362, 9762, 15677, 5618,
                                                                       5786, 11302, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17252, 3, 6122,
                                                                       6132, 11582, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17273, 3, 6132,
                                                                       6142, 11597, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17294, 3, 6142,
                                                                       6152, 11612, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17315, 3, 6152,
                                                                       6162, 11627, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17336, 3, 6162,
                                                                       6172, 11642, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17357, 3, 6172,
                                                                       6182, 11657, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17378, 3, 6182,
                                                                       6192, 11672, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17399, 3, 6192,
                                                                       6202, 11687, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17420, 3, 6202,
                                                                       6212, 11702, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17441, 0, 3,
                                                                       17252, 11582, 17273, 6232,
                                                                       6262, 11717, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17504, 0, 3,
                                                                       17273, 11597, 17294, 6262,
                                                                       6292, 11762, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17567, 0, 3,
                                                                       17294, 11612, 17315, 6292,
                                                                       6322, 11807, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17630, 0, 3,
                                                                       17315, 11627, 17336, 6322,
                                                                       6352, 11852, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17693, 0, 3,
                                                                       17336, 11642, 17357, 6352,
                                                                       6382, 11897, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17756, 0, 3,
                                                                       17357, 11657, 17378, 6382,
                                                                       6412, 11942, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17819, 0, 3,
                                                                       17378, 11672, 17399, 6412,
                                                                       6442, 11987, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17882, 0, 3,
                                                                       17399, 11687, 17420, 6442,
                                                                       6472, 12032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 17945, 0, 3,
                                                                       17441, 11717, 17504, 6532,
                                                                       6592, 12077, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 18071, 0, 3,
                                                                       17504, 11762, 17567, 6592,
                                                                       6652, 12167, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 18197, 0, 3,
                                                                       17567, 11807, 17630, 6652,
                                                                       6712, 12257, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 18323, 0, 3,
                                                                       17630, 11852, 17693, 6712,
                                                                       6772, 12347, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 18449, 0, 3,
                                                                       17693, 11897, 17756, 6772,
                                                                       6832, 12437, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 18575, 0, 3,
                                                                       17756, 11942, 17819, 6832,
                                                                       6892, 12527, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 18701, 0, 3,
                                                                       17819, 11987, 17882, 6892,
                                                                       6952, 12617, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 18827, 0, 3,
                                                                       17945, 12077, 18071, 7072,
                                                                       7172, 12707, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 19037, 0, 3,
                                                                       18071, 12167, 18197, 7172,
                                                                       7272, 12857, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 19247, 0, 3,
                                                                       18197, 12257, 18323, 7272,
                                                                       7372, 13007, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 19457, 0, 3,
                                                                       18323, 12347, 18449, 7372,
                                                                       7472, 13157, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 19667, 0, 3,
                                                                       18449, 12437, 18575, 7472,
                                                                       7572, 13307, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 19877, 0, 3,
                                                                       18575, 12527, 18701, 7572,
                                                                       7672, 13457, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 20087, 0, 3,
                                                                       18827, 12707, 19037, 7872,
                                                                       8022, 13607, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 20402, 0, 3,
                                                                       19037, 12857, 19247, 8022,
                                                                       8172, 13832, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 20717, 0, 3,
                                                                       19247, 13007, 19457, 8172,
                                                                       8322, 14057, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 21032, 0, 3,
                                                                       19457, 13157, 19667, 8322,
                                                                       8472, 14282, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 21347, 0, 3,
                                                                       19667, 13307, 19877, 8472,
                                                                       8622, 14507, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 21662, 0, 3,
                                                                       20087, 13607, 20402, 8922,
                                                                       9132, 14732, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 22103, 0, 3,
                                                                       20402, 13832, 20717, 9132,
                                                                       9342, 15047, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 22544, 0, 3,
                                                                       20717, 14057, 21032, 9342,
                                                                       9552, 15362, ncols, gamma,
                                                                       p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 22985, 0, 3,
                                                                       21032, 14282, 21347, 9552,
                                                                       9762, 15677, ncols, gamma,
                                                                       p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 23426, 0, 3,
                                                                       21662, 14732, 22103,
                                                                       10182, 10462, 15992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 24014, 0, 3,
                                                                       22103, 15047, 22544,
                                                                       10462, 10742, 16412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 24602, 0, 3,
                                                                       22544, 15362, 22985,
                                                                       10742, 11022, 16832,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25190, 3, 11582,
                                                                       11597, 17294, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25218, 3, 11597,
                                                                       11612, 17315, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25246, 3, 11612,
                                                                       11627, 17336, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25274, 3, 11627,
                                                                       11642, 17357, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25302, 3, 11642,
                                                                       11657, 17378, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25330, 3, 11657,
                                                                       11672, 17399, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 25358, 3, 11672,
                                                                       11687, 17420, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 25386, 0, 3,
                                                                       25190, 17294, 25218,
                                                                       11717, 11762, 17567,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 25470, 0, 3,
                                                                       25218, 17315, 25246,
                                                                       11762, 11807, 17630,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 25554, 0, 3,
                                                                       25246, 17336, 25274,
                                                                       11807, 11852, 17693,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 25638, 0, 3,
                                                                       25274, 17357, 25302,
                                                                       11852, 11897, 17756,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 25722, 0, 3,
                                                                       25302, 17378, 25330,
                                                                       11897, 11942, 17819,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 25806, 0, 3,
                                                                       25330, 17399, 25358,
                                                                       11942, 11987, 17882,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 25890, 0, 3,
                                                                       25386, 17567, 25470,
                                                                       12077, 12167, 18197,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 26058, 0, 3,
                                                                       25470, 17630, 25554,
                                                                       12167, 12257, 18323,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 26226, 0, 3,
                                                                       25554, 17693, 25638,
                                                                       12257, 12347, 18449,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 26394, 0, 3,
                                                                       25638, 17756, 25722,
                                                                       12347, 12437, 18575,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 26562, 0, 3,
                                                                       25722, 17819, 25806,
                                                                       12437, 12527, 18701,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 26730, 0, 3,
                                                                       25890, 18197, 26058,
                                                                       12707, 12857, 19247,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 27010, 0, 3,
                                                                       26058, 18323, 26226,
                                                                       12857, 13007, 19457,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 27290, 0, 3,
                                                                       26226, 18449, 26394,
                                                                       13007, 13157, 19667,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 27570, 0, 3,
                                                                       26394, 18575, 26562,
                                                                       13157, 13307, 19877,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 27850, 0, 3,
                                                                       26730, 19247, 27010,
                                                                       13607, 13832, 20717,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 28270, 0, 3,
                                                                       27010, 19457, 27290,
                                                                       13832, 14057, 21032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 28690, 0, 3,
                                                                       27290, 19667, 27570,
                                                                       14057, 14282, 21347,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 29110, 0, 3,
                                                                       27850, 20717, 28270,
                                                                       14732, 15047, 22544,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 29698, 0, 3,
                                                                       28270, 21032, 28690,
                                                                       15047, 15362, 22985,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 30286, 0, 3,
                                                                       29110, 22544, 29698,
                                                                       15992, 16412, 24602,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31070, 3, 17252,
                                                                       17273, 25190, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31106, 3, 17273,
                                                                       17294, 25218, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31142, 3, 17294,
                                                                       17315, 25246, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31178, 3, 17315,
                                                                       17336, 25274, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31214, 3, 17336,
                                                                       17357, 25302, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31250, 3, 17357,
                                                                       17378, 25330, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 31286, 3, 17378,
                                                                       17399, 25358, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 31322, 0, 3,
                                                                       31070, 25190, 31106,
                                                                       17441, 17504, 25386,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 31430, 0, 3,
                                                                       31106, 25218, 31142,
                                                                       17504, 17567, 25470,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 31538, 0, 3,
                                                                       31142, 25246, 31178,
                                                                       17567, 17630, 25554,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 31646, 0, 3,
                                                                       31178, 25274, 31214,
                                                                       17630, 17693, 25638,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 31754, 0, 3,
                                                                       31214, 25302, 31250,
                                                                       17693, 17756, 25722,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 31862, 0, 3,
                                                                       31250, 25330, 31286,
                                                                       17756, 17819, 25806,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 31970, 0, 3,
                                                                       31322, 25386, 31430,
                                                                       17945, 18071, 25890,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 32186, 0, 3,
                                                                       31430, 25470, 31538,
                                                                       18071, 18197, 26058,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 32402, 0, 3,
                                                                       31538, 25554, 31646,
                                                                       18197, 18323, 26226,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 32618, 0, 3,
                                                                       31646, 25638, 31754,
                                                                       18323, 18449, 26394,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 32834, 0, 3,
                                                                       31754, 25722, 31862,
                                                                       18449, 18575, 26562,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 33050, 0, 3,
                                                                       31970, 25890, 32186,
                                                                       18827, 19037, 26730,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 33410, 0, 3,
                                                                       32186, 26058, 32402,
                                                                       19037, 19247, 27010,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 33770, 0, 3,
                                                                       32402, 26226, 32618,
                                                                       19247, 19457, 27290,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 34130, 0, 3,
                                                                       32618, 26394, 32834,
                                                                       19457, 19667, 27570,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 34490, 0, 3,
                                                                       33050, 26730, 33410,
                                                                       20087, 20402, 27850,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 35030, 0, 3,
                                                                       33410, 27010, 33770,
                                                                       20402, 20717, 28270,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 35570, 0, 3,
                                                                       33770, 27290, 34130,
                                                                       20717, 21032, 28690,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 36110, 0, 3,
                                                                       34490, 27850, 35030,
                                                                       21662, 22103, 29110,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 36866, 0, 3,
                                                                       35030, 28270, 35570,
                                                                       22103, 22544, 29698,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 37622, 0, 3,
                                                                       36110, 29110, 36866,
                                                                       23426, 24014, 30286,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 38630, 36110, 756, ncols);

                    simdfunc::contract_primitives(buffer, 39701, 37622, 1008, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 39386, 38630, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 40709, 39701, 28, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 41129, 39386, 40709, 15, nmax);

        simdtrf::transform_h_inner(buffer, 42074, 41129, 3, 15, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 42074, 165, nmax);
    }

    for (size_t m = 0; m < 495; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
