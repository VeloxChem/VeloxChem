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


#include "SimdThreeCenterElectronRepulsionRecFHH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDH.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferFH.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_fhh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_fhh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 49349, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 847 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 49349, 35630, 3665, dimensions);

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

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 722, 0, 3, 358,
                                                                       379, 526, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 758, 0, 3, 379,
                                                                       400, 554, 582, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 794, 0, 3, 400,
                                                                       421, 582, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 830, 0, 3, 421,
                                                                       442, 610, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 866, 0, 3, 442,
                                                                       463, 638, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 902, 0, 3, 463,
                                                                       484, 666, 694, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 938, 0, 3, 526,
                                                                       554, 722, 758, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 983, 0, 3, 554,
                                                                       582, 758, 794, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1028, 0, 3, 582,
                                                                       610, 794, 830, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1073, 0, 3, 610,
                                                                       638, 830, 866, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1118, 0, 3, 638,
                                                                       666, 866, 902, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1163, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1166, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1169, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1172, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1175, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1178, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1181, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1184, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1187, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1190, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1193, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1196, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1199, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1202, 3, 10, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1211, 3, 11, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1220, 3, 12, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1229, 3, 13, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1238, 3, 14, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1247, 3, 15, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1256, 3, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1265, 3, 17, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1274, 3, 18, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1283, 3, 19, 54,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1292, 3, 21, 57,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1310, 3, 24, 63,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1328, 3, 27, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1346, 3, 30, 75,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1364, 3, 33, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1382, 3, 36, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1400, 3, 39, 93,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1418, 3, 42, 99,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1436, 3, 45, 105,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1454, 3, 48, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1472, 3, 51, 117,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1490, 3, 57, 123,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1520, 3, 63, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1550, 3, 69, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1580, 3, 75, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1610, 3, 81, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1640, 3, 87, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1670, 3, 93, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1700, 3, 99, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1730, 3, 105, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1760, 3, 111, 213,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1790, 3, 123, 223,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1835, 3, 133, 238,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1880, 3, 143, 253,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1925, 3, 153, 268,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1970, 3, 163, 283,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2015, 3, 173, 298,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2060, 3, 183, 313,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2105, 3, 193, 328,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2150, 3, 203, 343,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2195, 3, 223, 358,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2258, 3, 238, 379,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2321, 3, 253, 400,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2384, 3, 268, 421,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2447, 3, 283, 442,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2510, 3, 298, 463,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2573, 3, 313, 484,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2636, 3, 328, 505,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2699, 3, 358, 526,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2783, 3, 379, 554,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2867, 3, 400, 582,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2951, 3, 421, 610,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3035, 3, 442, 638,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3119, 3, 463, 666,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3203, 3, 484, 694,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3287, 3, 526, 722,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3395, 3, 554, 758,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3503, 3, 582, 794,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3611, 3, 610, 830,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3719, 3, 638, 866,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3827, 3, 666, 902,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3935, 3, 722, 938,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4070, 3, 758, 983,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4205, 3, 794,
                                                                       1028, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4340, 3, 830,
                                                                       1073, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4475, 3, 866,
                                                                       1118, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4610, 3, 8, 9,
                                                                       1169, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4616, 3, 9, 10,
                                                                       1172, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4622, 3, 10, 11,
                                                                       1175, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4628, 3, 11, 12,
                                                                       1178, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4634, 3, 12, 13,
                                                                       1181, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4640, 3, 13, 14,
                                                                       1184, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4646, 3, 14, 15,
                                                                       1187, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4652, 3, 15, 16,
                                                                       1190, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4658, 3, 16, 17,
                                                                       1193, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4664, 3, 17, 18,
                                                                       1196, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4670, 3, 18, 19,
                                                                       1199, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4676, 0, 3, 4610,
                                                                       1169, 4616, 1202, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4694, 0, 3, 4616,
                                                                       1172, 4622, 1211, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4712, 0, 3, 4622,
                                                                       1175, 4628, 1220, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4730, 0, 3, 4628,
                                                                       1178, 4634, 1229, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4748, 0, 3, 4634,
                                                                       1181, 4640, 1238, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4766, 0, 3, 4640,
                                                                       1184, 4646, 1247, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4784, 0, 3, 4646,
                                                                       1187, 4652, 1256, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4802, 0, 3, 4652,
                                                                       1190, 4658, 1265, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4820, 0, 3, 4658,
                                                                       1193, 4664, 1274, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4838, 0, 3, 4664,
                                                                       1196, 4670, 1283, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4856, 0, 3, 4676,
                                                                       1202, 4694, 57, 63, 1328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4892, 0, 3, 4694,
                                                                       1211, 4712, 63, 69, 1346,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4928, 0, 3, 4712,
                                                                       1220, 4730, 69, 75, 1364,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4964, 0, 3, 4730,
                                                                       1229, 4748, 75, 81, 1382,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5000, 0, 3, 4748,
                                                                       1238, 4766, 81, 87, 1400,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5036, 0, 3, 4766,
                                                                       1247, 4784, 87, 93, 1418,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5072, 0, 3, 4784,
                                                                       1256, 4802, 93, 99, 1436,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5108, 0, 3, 4802,
                                                                       1265, 4820, 99, 105, 1454,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5144, 0, 3, 4820,
                                                                       1274, 4838, 105, 111,
                                                                       1472, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5180, 0, 3, 4856,
                                                                       1328, 4892, 123, 133,
                                                                       1550, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5240, 0, 3, 4892,
                                                                       1346, 4928, 133, 143,
                                                                       1580, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5300, 0, 3, 4928,
                                                                       1364, 4964, 143, 153,
                                                                       1610, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5360, 0, 3, 4964,
                                                                       1382, 5000, 153, 163,
                                                                       1640, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5420, 0, 3, 5000,
                                                                       1400, 5036, 163, 173,
                                                                       1670, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5480, 0, 3, 5036,
                                                                       1418, 5072, 173, 183,
                                                                       1700, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5540, 0, 3, 5072,
                                                                       1436, 5108, 183, 193,
                                                                       1730, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5600, 0, 3, 5108,
                                                                       1454, 5144, 193, 203,
                                                                       1760, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5660, 0, 3, 5180,
                                                                       1550, 5240, 223, 238,
                                                                       1880, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5750, 0, 3, 5240,
                                                                       1580, 5300, 238, 253,
                                                                       1925, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5840, 0, 3, 5300,
                                                                       1610, 5360, 253, 268,
                                                                       1970, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5930, 0, 3, 5360,
                                                                       1640, 5420, 268, 283,
                                                                       2015, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6020, 0, 3, 5420,
                                                                       1670, 5480, 283, 298,
                                                                       2060, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6110, 0, 3, 5480,
                                                                       1700, 5540, 298, 313,
                                                                       2105, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6200, 0, 3, 5540,
                                                                       1730, 5600, 313, 328,
                                                                       2150, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6290, 0, 3, 5660,
                                                                       1880, 5750, 358, 379,
                                                                       2321, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6416, 0, 3, 5750,
                                                                       1925, 5840, 379, 400,
                                                                       2384, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6542, 0, 3, 5840,
                                                                       1970, 5930, 400, 421,
                                                                       2447, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6668, 0, 3, 5930,
                                                                       2015, 6020, 421, 442,
                                                                       2510, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6794, 0, 3, 6020,
                                                                       2060, 6110, 442, 463,
                                                                       2573, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6920, 0, 3, 6110,
                                                                       2105, 6200, 463, 484,
                                                                       2636, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7046, 0, 3, 6290,
                                                                       2321, 6416, 526, 554,
                                                                       2867, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7214, 0, 3, 6416,
                                                                       2384, 6542, 554, 582,
                                                                       2951, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7382, 0, 3, 6542,
                                                                       2447, 6668, 582, 610,
                                                                       3035, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7550, 0, 3, 6668,
                                                                       2510, 6794, 610, 638,
                                                                       3119, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7718, 0, 3, 6794,
                                                                       2573, 6920, 638, 666,
                                                                       3203, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 7886, 0, 3, 7046,
                                                                       2867, 7214, 722, 758,
                                                                       3503, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 8102, 0, 3, 7214,
                                                                       2951, 7382, 758, 794,
                                                                       3611, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 8318, 0, 3, 7382,
                                                                       3035, 7550, 794, 830,
                                                                       3719, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 8534, 0, 3, 7550,
                                                                       3119, 7718, 830, 866,
                                                                       3827, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 8750, 0, 3, 7886,
                                                                       3503, 8102, 938, 983,
                                                                       4205, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 9020, 0, 3, 8102,
                                                                       3611, 8318, 983, 1028,
                                                                       4340, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 9290, 0, 3, 8318,
                                                                       3719, 8534, 1028, 1073,
                                                                       4475, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9560, 3, 1163,
                                                                       1166, 4610, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9570, 3, 1166,
                                                                       1169, 4616, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9580, 3, 1169,
                                                                       1172, 4622, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9590, 3, 1172,
                                                                       1175, 4628, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9600, 3, 1175,
                                                                       1178, 4634, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9610, 3, 1178,
                                                                       1181, 4640, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9620, 3, 1181,
                                                                       1184, 4646, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9630, 3, 1184,
                                                                       1187, 4652, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9640, 3, 1187,
                                                                       1190, 4658, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9650, 3, 1190,
                                                                       1193, 4664, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9660, 3, 1193,
                                                                       1196, 4670, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9670, 0, 3, 9560,
                                                                       4610, 9570, 4676, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9700, 0, 3, 9570,
                                                                       4616, 9580, 4694, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9730, 0, 3, 9580,
                                                                       4622, 9590, 4712, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9760, 0, 3, 9590,
                                                                       4628, 9600, 4730, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9790, 0, 3, 9600,
                                                                       4634, 9610, 4748, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9820, 0, 3, 9610,
                                                                       4640, 9620, 4766, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9850, 0, 3, 9620,
                                                                       4646, 9630, 4784, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9880, 0, 3, 9630,
                                                                       4652, 9640, 4802, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9910, 0, 3, 9640,
                                                                       4658, 9650, 4820, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9940, 0, 3, 9650,
                                                                       4664, 9660, 4838, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9970, 0, 3, 9670,
                                                                       4676, 9700, 1292, 1310,
                                                                       4856, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10030, 0, 3, 9700,
                                                                       4694, 9730, 1310, 1328,
                                                                       4892, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10090, 0, 3, 9730,
                                                                       4712, 9760, 1328, 1346,
                                                                       4928, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10150, 0, 3, 9760,
                                                                       4730, 9790, 1346, 1364,
                                                                       4964, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10210, 0, 3, 9790,
                                                                       4748, 9820, 1364, 1382,
                                                                       5000, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10270, 0, 3, 9820,
                                                                       4766, 9850, 1382, 1400,
                                                                       5036, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10330, 0, 3, 9850,
                                                                       4784, 9880, 1400, 1418,
                                                                       5072, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10390, 0, 3, 9880,
                                                                       4802, 9910, 1418, 1436,
                                                                       5108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10450, 0, 3, 9910,
                                                                       4820, 9940, 1436, 1454,
                                                                       5144, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10510, 0, 3, 9970,
                                                                       4856, 10030, 1490, 1520,
                                                                       5180, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10610, 0, 3,
                                                                       10030, 4892, 10090, 1520,
                                                                       1550, 5240, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10710, 0, 3,
                                                                       10090, 4928, 10150, 1550,
                                                                       1580, 5300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10810, 0, 3,
                                                                       10150, 4964, 10210, 1580,
                                                                       1610, 5360, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10910, 0, 3,
                                                                       10210, 5000, 10270, 1610,
                                                                       1640, 5420, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11010, 0, 3,
                                                                       10270, 5036, 10330, 1640,
                                                                       1670, 5480, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11110, 0, 3,
                                                                       10330, 5072, 10390, 1670,
                                                                       1700, 5540, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11210, 0, 3,
                                                                       10390, 5108, 10450, 1700,
                                                                       1730, 5600, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11310, 0, 3,
                                                                       10510, 5180, 10610, 1790,
                                                                       1835, 5660, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11460, 0, 3,
                                                                       10610, 5240, 10710, 1835,
                                                                       1880, 5750, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11610, 0, 3,
                                                                       10710, 5300, 10810, 1880,
                                                                       1925, 5840, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11760, 0, 3,
                                                                       10810, 5360, 10910, 1925,
                                                                       1970, 5930, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11910, 0, 3,
                                                                       10910, 5420, 11010, 1970,
                                                                       2015, 6020, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12060, 0, 3,
                                                                       11010, 5480, 11110, 2015,
                                                                       2060, 6110, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12210, 0, 3,
                                                                       11110, 5540, 11210, 2060,
                                                                       2105, 6200, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12360, 0, 3,
                                                                       11310, 5660, 11460, 2195,
                                                                       2258, 6290, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12570, 0, 3,
                                                                       11460, 5750, 11610, 2258,
                                                                       2321, 6416, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12780, 0, 3,
                                                                       11610, 5840, 11760, 2321,
                                                                       2384, 6542, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12990, 0, 3,
                                                                       11760, 5930, 11910, 2384,
                                                                       2447, 6668, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 13200, 0, 3,
                                                                       11910, 6020, 12060, 2447,
                                                                       2510, 6794, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 13410, 0, 3,
                                                                       12060, 6110, 12210, 2510,
                                                                       2573, 6920, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 13620, 0, 3,
                                                                       12360, 6290, 12570, 2699,
                                                                       2783, 7046, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 13900, 0, 3,
                                                                       12570, 6416, 12780, 2783,
                                                                       2867, 7214, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 14180, 0, 3,
                                                                       12780, 6542, 12990, 2867,
                                                                       2951, 7382, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 14460, 0, 3,
                                                                       12990, 6668, 13200, 2951,
                                                                       3035, 7550, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 14740, 0, 3,
                                                                       13200, 6794, 13410, 3035,
                                                                       3119, 7718, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 15020, 0, 3,
                                                                       13620, 7046, 13900, 3287,
                                                                       3395, 7886, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 15380, 0, 3,
                                                                       13900, 7214, 14180, 3395,
                                                                       3503, 8102, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 15740, 0, 3,
                                                                       14180, 7382, 14460, 3503,
                                                                       3611, 8318, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 16100, 0, 3,
                                                                       14460, 7550, 14740, 3611,
                                                                       3719, 8534, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 16460, 0, 3,
                                                                       15020, 7886, 15380, 3935,
                                                                       4070, 8750, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 16910, 0, 3,
                                                                       15380, 8102, 15740, 4070,
                                                                       4205, 9020, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 17360, 0, 3,
                                                                       15740, 8318, 16100, 4205,
                                                                       4340, 9290, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17810, 3, 4610,
                                                                       4616, 9580, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17825, 3, 4616,
                                                                       4622, 9590, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17840, 3, 4622,
                                                                       4628, 9600, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17855, 3, 4628,
                                                                       4634, 9610, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17870, 3, 4634,
                                                                       4640, 9620, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17885, 3, 4640,
                                                                       4646, 9630, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17900, 3, 4646,
                                                                       4652, 9640, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17915, 3, 4652,
                                                                       4658, 9650, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17930, 3, 4658,
                                                                       4664, 9660, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17945, 0, 3,
                                                                       17810, 9580, 17825, 4676,
                                                                       4694, 9730, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17990, 0, 3,
                                                                       17825, 9590, 17840, 4694,
                                                                       4712, 9760, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18035, 0, 3,
                                                                       17840, 9600, 17855, 4712,
                                                                       4730, 9790, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18080, 0, 3,
                                                                       17855, 9610, 17870, 4730,
                                                                       4748, 9820, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18125, 0, 3,
                                                                       17870, 9620, 17885, 4748,
                                                                       4766, 9850, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18170, 0, 3,
                                                                       17885, 9630, 17900, 4766,
                                                                       4784, 9880, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18215, 0, 3,
                                                                       17900, 9640, 17915, 4784,
                                                                       4802, 9910, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18260, 0, 3,
                                                                       17915, 9650, 17930, 4802,
                                                                       4820, 9940, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18305, 0, 3,
                                                                       17945, 9730, 17990, 4856,
                                                                       4892, 10090, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18395, 0, 3,
                                                                       17990, 9760, 18035, 4892,
                                                                       4928, 10150, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18485, 0, 3,
                                                                       18035, 9790, 18080, 4928,
                                                                       4964, 10210, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18575, 0, 3,
                                                                       18080, 9820, 18125, 4964,
                                                                       5000, 10270, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18665, 0, 3,
                                                                       18125, 9850, 18170, 5000,
                                                                       5036, 10330, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18755, 0, 3,
                                                                       18170, 9880, 18215, 5036,
                                                                       5072, 10390, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18845, 0, 3,
                                                                       18215, 9910, 18260, 5072,
                                                                       5108, 10450, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 18935, 0, 3,
                                                                       18305, 10090, 18395, 5180,
                                                                       5240, 10710, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19085, 0, 3,
                                                                       18395, 10150, 18485, 5240,
                                                                       5300, 10810, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19235, 0, 3,
                                                                       18485, 10210, 18575, 5300,
                                                                       5360, 10910, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19385, 0, 3,
                                                                       18575, 10270, 18665, 5360,
                                                                       5420, 11010, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19535, 0, 3,
                                                                       18665, 10330, 18755, 5420,
                                                                       5480, 11110, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19685, 0, 3,
                                                                       18755, 10390, 18845, 5480,
                                                                       5540, 11210, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 19835, 0, 3,
                                                                       18935, 10710, 19085, 5660,
                                                                       5750, 11610, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 20060, 0, 3,
                                                                       19085, 10810, 19235, 5750,
                                                                       5840, 11760, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 20285, 0, 3,
                                                                       19235, 10910, 19385, 5840,
                                                                       5930, 11910, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 20510, 0, 3,
                                                                       19385, 11010, 19535, 5930,
                                                                       6020, 12060, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 20735, 0, 3,
                                                                       19535, 11110, 19685, 6020,
                                                                       6110, 12210, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 20960, 0, 3,
                                                                       19835, 11610, 20060, 6290,
                                                                       6416, 12780, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 21275, 0, 3,
                                                                       20060, 11760, 20285, 6416,
                                                                       6542, 12990, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 21590, 0, 3,
                                                                       20285, 11910, 20510, 6542,
                                                                       6668, 13200, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 21905, 0, 3,
                                                                       20510, 12060, 20735, 6668,
                                                                       6794, 13410, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 22220, 0, 3,
                                                                       20960, 12780, 21275, 7046,
                                                                       7214, 14180, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 22640, 0, 3,
                                                                       21275, 12990, 21590, 7214,
                                                                       7382, 14460, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 23060, 0, 3,
                                                                       21590, 13200, 21905, 7382,
                                                                       7550, 14740, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 23480, 0, 3,
                                                                       22220, 14180, 22640, 7886,
                                                                       8102, 15740, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 24020, 0, 3,
                                                                       22640, 14460, 23060, 8102,
                                                                       8318, 16100, ncols, gamma,
                                                                       p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 24560, 0, 3,
                                                                       23480, 15740, 24020, 8750,
                                                                       9020, 17360, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25235, 3, 9560,
                                                                       9570, 17810, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25256, 3, 9570,
                                                                       9580, 17825, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25277, 3, 9580,
                                                                       9590, 17840, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25298, 3, 9590,
                                                                       9600, 17855, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25319, 3, 9600,
                                                                       9610, 17870, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25340, 3, 9610,
                                                                       9620, 17885, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25361, 3, 9620,
                                                                       9630, 17900, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25382, 3, 9630,
                                                                       9640, 17915, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25403, 3, 9640,
                                                                       9650, 17930, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 25424, 0, 3,
                                                                       25235, 17810, 25256, 9670,
                                                                       9700, 17945, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 25487, 0, 3,
                                                                       25256, 17825, 25277, 9700,
                                                                       9730, 17990, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 25550, 0, 3,
                                                                       25277, 17840, 25298, 9730,
                                                                       9760, 18035, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 25613, 0, 3,
                                                                       25298, 17855, 25319, 9760,
                                                                       9790, 18080, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 25676, 0, 3,
                                                                       25319, 17870, 25340, 9790,
                                                                       9820, 18125, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 25739, 0, 3,
                                                                       25340, 17885, 25361, 9820,
                                                                       9850, 18170, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 25802, 0, 3,
                                                                       25361, 17900, 25382, 9850,
                                                                       9880, 18215, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 25865, 0, 3,
                                                                       25382, 17915, 25403, 9880,
                                                                       9910, 18260, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 25928, 0, 3,
                                                                       25424, 17945, 25487, 9970,
                                                                       10030, 18305, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26054, 0, 3,
                                                                       25487, 17990, 25550,
                                                                       10030, 10090, 18395,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26180, 0, 3,
                                                                       25550, 18035, 25613,
                                                                       10090, 10150, 18485,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26306, 0, 3,
                                                                       25613, 18080, 25676,
                                                                       10150, 10210, 18575,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26432, 0, 3,
                                                                       25676, 18125, 25739,
                                                                       10210, 10270, 18665,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26558, 0, 3,
                                                                       25739, 18170, 25802,
                                                                       10270, 10330, 18755,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26684, 0, 3,
                                                                       25802, 18215, 25865,
                                                                       10330, 10390, 18845,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 26810, 0, 3,
                                                                       25928, 18305, 26054,
                                                                       10510, 10610, 18935,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 27020, 0, 3,
                                                                       26054, 18395, 26180,
                                                                       10610, 10710, 19085,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 27230, 0, 3,
                                                                       26180, 18485, 26306,
                                                                       10710, 10810, 19235,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 27440, 0, 3,
                                                                       26306, 18575, 26432,
                                                                       10810, 10910, 19385,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 27650, 0, 3,
                                                                       26432, 18665, 26558,
                                                                       10910, 11010, 19535,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 27860, 0, 3,
                                                                       26558, 18755, 26684,
                                                                       11010, 11110, 19685,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 28070, 0, 3,
                                                                       26810, 18935, 27020,
                                                                       11310, 11460, 19835,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 28385, 0, 3,
                                                                       27020, 19085, 27230,
                                                                       11460, 11610, 20060,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 28700, 0, 3,
                                                                       27230, 19235, 27440,
                                                                       11610, 11760, 20285,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 29015, 0, 3,
                                                                       27440, 19385, 27650,
                                                                       11760, 11910, 20510,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 29330, 0, 3,
                                                                       27650, 19535, 27860,
                                                                       11910, 12060, 20735,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 29645, 0, 3,
                                                                       28070, 19835, 28385,
                                                                       12360, 12570, 20960,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 30086, 0, 3,
                                                                       28385, 20060, 28700,
                                                                       12570, 12780, 21275,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 30527, 0, 3,
                                                                       28700, 20285, 29015,
                                                                       12780, 12990, 21590,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 30968, 0, 3,
                                                                       29015, 20510, 29330,
                                                                       12990, 13200, 21905,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 31409, 0, 3,
                                                                       29645, 20960, 30086,
                                                                       13620, 13900, 22220,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 31997, 0, 3,
                                                                       30086, 21275, 30527,
                                                                       13900, 14180, 22640,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 32585, 0, 3,
                                                                       30527, 21590, 30968,
                                                                       14180, 14460, 23060,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 33173, 0, 3,
                                                                       31409, 22220, 31997,
                                                                       15020, 15380, 23480,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 33929, 0, 3,
                                                                       31997, 22640, 32585,
                                                                       15380, 15740, 24020,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 34685, 0, 3,
                                                                       33173, 23480, 33929,
                                                                       16460, 16910, 24560,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 35630, 29645, 441, ncols);

                    simdfunc::contract_primitives(buffer, 36302, 31409, 588, ncols);

                    simdfunc::contract_primitives(buffer, 37198, 33173, 756, ncols);

                    simdfunc::contract_primitives(buffer, 38350, 34685, 945, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 36071, 35630, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 36890, 36302, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 37954, 37198, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 39295, 38350, 45, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 39790, 36071, 36890, 11, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 40483, 36890, 37954, 11, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 41407, 37954, 39295, 11, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 42595, 39790, 40483, 11, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 43981, 40483, 41407, 11, nmax);

        simdtrf::compute_hrr_fh(buffer, coordinates, 45829, 42595, 43981, 11, nmax);

        simdtrf::transform_h_inner(buffer, 48139, 45829, 10, 11, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 48139, 121, nmax);
    }

    for (size_t m = 0; m < 847; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
