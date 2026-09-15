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


#include "SimdThreeCenterElectronRepulsionRecSHH.hpp"

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
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_shh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_shh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 10406, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 121 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 10406, 9734, 441, dimensions);

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
                                                        5, 6, 7, 8, 9, 10}, ncols, fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 18, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 45, 0, 3, 8, 9,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 51, 0, 3, 9, 10,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 57, 0, 3, 10, 11,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 63, 0, 3, 11, 12,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 69, 0, 3, 12, 13,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 75, 0, 3, 13, 14,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 81, 0, 3, 14, 15,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 87, 0, 3, 15, 16,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 93, 0, 3, 18, 21,
                                                                       45, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 103, 0, 3, 21, 24,
                                                                       51, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 113, 0, 3, 24, 27,
                                                                       57, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 123, 0, 3, 27, 30,
                                                                       63, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 133, 0, 3, 30, 33,
                                                                       69, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 143, 0, 3, 33, 36,
                                                                       75, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 153, 0, 3, 36, 39,
                                                                       81, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 163, 0, 3, 45, 51,
                                                                       93, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 178, 0, 3, 51, 57,
                                                                       103, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 193, 0, 3, 57, 63,
                                                                       113, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 208, 0, 3, 63, 69,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 223, 0, 3, 69, 75,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 75, 81,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 253, 0, 3, 93,
                                                                       103, 163, 178, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 274, 0, 3, 103,
                                                                       113, 178, 193, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 295, 0, 3, 113,
                                                                       123, 193, 208, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 316, 0, 3, 123,
                                                                       133, 208, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 337, 0, 3, 133,
                                                                       143, 223, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 358, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 361, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 364, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 367, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 370, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 373, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 376, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 379, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 382, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 385, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 388, 3, 10, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 397, 3, 11, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 406, 3, 12, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 415, 3, 13, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 424, 3, 14, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 433, 3, 15, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 442, 3, 16, 42,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 451, 3, 18, 45,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 469, 3, 21, 51,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 487, 3, 24, 57,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 505, 3, 27, 63,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 523, 3, 30, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 541, 3, 33, 75,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 559, 3, 36, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 577, 3, 39, 87,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 595, 3, 45, 93,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 625, 3, 51, 103,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 655, 3, 57, 113,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 685, 3, 63, 123,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 715, 3, 69, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 745, 3, 75, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 775, 3, 81, 153,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 805, 3, 93, 163,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 850, 3, 103, 178,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 895, 3, 113, 193,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 940, 3, 123, 208,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 985, 3, 133, 223,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1030, 3, 143, 238,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1075, 3, 163, 253,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1138, 3, 178, 274,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1201, 3, 193, 295,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1264, 3, 208, 316,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1327, 3, 223, 337,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1390, 3, 8, 9,
                                                                       364, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1396, 3, 9, 10,
                                                                       367, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1402, 3, 10, 11,
                                                                       370, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1408, 3, 11, 12,
                                                                       373, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1414, 3, 12, 13,
                                                                       376, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1420, 3, 13, 14,
                                                                       379, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1426, 3, 14, 15,
                                                                       382, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1432, 3, 15, 16,
                                                                       385, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1438, 0, 3, 1390,
                                                                       364, 1396, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1456, 0, 3, 1396,
                                                                       367, 1402, 397, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1474, 0, 3, 1402,
                                                                       370, 1408, 406, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 1408,
                                                                       373, 1414, 415, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1510, 0, 3, 1414,
                                                                       376, 1420, 424, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1528, 0, 3, 1420,
                                                                       379, 1426, 433, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1546, 0, 3, 1426,
                                                                       382, 1432, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1564, 0, 3, 1438,
                                                                       388, 1456, 45, 51, 487,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1600, 0, 3, 1456,
                                                                       397, 1474, 51, 57, 505,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1636, 0, 3, 1474,
                                                                       406, 1492, 57, 63, 523,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1672, 0, 3, 1492,
                                                                       415, 1510, 63, 69, 541,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1708, 0, 3, 1510,
                                                                       424, 1528, 69, 75, 559,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 1528,
                                                                       433, 1546, 75, 81, 577,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1780, 0, 3, 1564,
                                                                       487, 1600, 93, 103, 655,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1840, 0, 3, 1600,
                                                                       505, 1636, 103, 113, 685,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1900, 0, 3, 1636,
                                                                       523, 1672, 113, 123, 715,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1960, 0, 3, 1672,
                                                                       541, 1708, 123, 133, 745,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2020, 0, 3, 1708,
                                                                       559, 1744, 133, 143, 775,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2080, 0, 3, 1780,
                                                                       655, 1840, 163, 178, 895,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2170, 0, 3, 1840,
                                                                       685, 1900, 178, 193, 940,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2260, 0, 3, 1900,
                                                                       715, 1960, 193, 208, 985,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2350, 0, 3, 1960,
                                                                       745, 2020, 208, 223, 1030,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2440, 0, 3, 2080,
                                                                       895, 2170, 253, 274, 1201,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2566, 0, 3, 2170,
                                                                       940, 2260, 274, 295, 1264,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2692, 0, 3, 2260,
                                                                       985, 2350, 295, 316, 1327,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2818, 3, 358, 361,
                                                                       1390, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2828, 3, 361, 364,
                                                                       1396, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2838, 3, 364, 367,
                                                                       1402, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2848, 3, 367, 370,
                                                                       1408, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2858, 3, 370, 373,
                                                                       1414, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2868, 3, 373, 376,
                                                                       1420, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2878, 3, 376, 379,
                                                                       1426, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2888, 3, 379, 382,
                                                                       1432, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2898, 0, 3, 2818,
                                                                       1390, 2828, 1438, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2928, 0, 3, 2828,
                                                                       1396, 2838, 1456, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2958, 0, 3, 2838,
                                                                       1402, 2848, 1474, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2988, 0, 3, 2848,
                                                                       1408, 2858, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3018, 0, 3, 2858,
                                                                       1414, 2868, 1510, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3048, 0, 3, 2868,
                                                                       1420, 2878, 1528, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3078, 0, 3, 2878,
                                                                       1426, 2888, 1546, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3108, 0, 3, 2898,
                                                                       1438, 2928, 451, 469,
                                                                       1564, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3168, 0, 3, 2928,
                                                                       1456, 2958, 469, 487,
                                                                       1600, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3228, 0, 3, 2958,
                                                                       1474, 2988, 487, 505,
                                                                       1636, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3288, 0, 3, 2988,
                                                                       1492, 3018, 505, 523,
                                                                       1672, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3348, 0, 3, 3018,
                                                                       1510, 3048, 523, 541,
                                                                       1708, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3408, 0, 3, 3048,
                                                                       1528, 3078, 541, 559,
                                                                       1744, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3468, 0, 3, 3108,
                                                                       1564, 3168, 595, 625,
                                                                       1780, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3568, 0, 3, 3168,
                                                                       1600, 3228, 625, 655,
                                                                       1840, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3668, 0, 3, 3228,
                                                                       1636, 3288, 655, 685,
                                                                       1900, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3768, 0, 3, 3288,
                                                                       1672, 3348, 685, 715,
                                                                       1960, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 3868, 0, 3, 3348,
                                                                       1708, 3408, 715, 745,
                                                                       2020, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 3468,
                                                                       1780, 3568, 805, 850,
                                                                       2080, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 4118, 0, 3, 3568,
                                                                       1840, 3668, 850, 895,
                                                                       2170, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 4268, 0, 3, 3668,
                                                                       1900, 3768, 895, 940,
                                                                       2260, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 4418, 0, 3, 3768,
                                                                       1960, 3868, 940, 985,
                                                                       2350, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 4568, 0, 3, 3968,
                                                                       2080, 4118, 1075, 1138,
                                                                       2440, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 4778, 0, 3, 4118,
                                                                       2170, 4268, 1138, 1201,
                                                                       2566, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 4988, 0, 3, 4268,
                                                                       2260, 4418, 1201, 1264,
                                                                       2692, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5198, 3, 1390,
                                                                       1396, 2838, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5213, 3, 1396,
                                                                       1402, 2848, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5228, 3, 1402,
                                                                       1408, 2858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5243, 3, 1408,
                                                                       1414, 2868, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5258, 3, 1414,
                                                                       1420, 2878, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5273, 3, 1420,
                                                                       1426, 2888, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 5198,
                                                                       2838, 5213, 1438, 1456,
                                                                       2958, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5333, 0, 3, 5213,
                                                                       2848, 5228, 1456, 1474,
                                                                       2988, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5378, 0, 3, 5228,
                                                                       2858, 5243, 1474, 1492,
                                                                       3018, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5423, 0, 3, 5243,
                                                                       2868, 5258, 1492, 1510,
                                                                       3048, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 5468, 0, 3, 5258,
                                                                       2878, 5273, 1510, 1528,
                                                                       3078, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5513, 0, 3, 5288,
                                                                       2958, 5333, 1564, 1600,
                                                                       3228, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5603, 0, 3, 5333,
                                                                       2988, 5378, 1600, 1636,
                                                                       3288, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5693, 0, 3, 5378,
                                                                       3018, 5423, 1636, 1672,
                                                                       3348, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5783, 0, 3, 5423,
                                                                       3048, 5468, 1672, 1708,
                                                                       3408, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 5873, 0, 3, 5513,
                                                                       3228, 5603, 1780, 1840,
                                                                       3668, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 6023, 0, 3, 5603,
                                                                       3288, 5693, 1840, 1900,
                                                                       3768, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 6173, 0, 3, 5693,
                                                                       3348, 5783, 1900, 1960,
                                                                       3868, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 6323, 0, 3, 5873,
                                                                       3668, 6023, 2080, 2170,
                                                                       4268, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 6548, 0, 3, 6023,
                                                                       3768, 6173, 2170, 2260,
                                                                       4418, ncols, gamma, p,
                                                                       q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 6773, 0, 3, 6323,
                                                                       4268, 6548, 2440, 2566,
                                                                       4988, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7088, 3, 2818,
                                                                       2828, 5198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7109, 3, 2828,
                                                                       2838, 5213, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7130, 3, 2838,
                                                                       2848, 5228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7151, 3, 2848,
                                                                       2858, 5243, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7172, 3, 2858,
                                                                       2868, 5258, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7193, 3, 2868,
                                                                       2878, 5273, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7214, 0, 3, 7088,
                                                                       5198, 7109, 2898, 2928,
                                                                       5288, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7277, 0, 3, 7109,
                                                                       5213, 7130, 2928, 2958,
                                                                       5333, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7340, 0, 3, 7130,
                                                                       5228, 7151, 2958, 2988,
                                                                       5378, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7403, 0, 3, 7151,
                                                                       5243, 7172, 2988, 3018,
                                                                       5423, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 7466, 0, 3, 7172,
                                                                       5258, 7193, 3018, 3048,
                                                                       5468, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 7529, 0, 3, 7214,
                                                                       5288, 7277, 3108, 3168,
                                                                       5513, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 7655, 0, 3, 7277,
                                                                       5333, 7340, 3168, 3228,
                                                                       5603, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 7781, 0, 3, 7340,
                                                                       5378, 7403, 3228, 3288,
                                                                       5693, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 7907, 0, 3, 7403,
                                                                       5423, 7466, 3288, 3348,
                                                                       5783, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 8033, 0, 3, 7529,
                                                                       5513, 7655, 3468, 3568,
                                                                       5873, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 8243, 0, 3, 7655,
                                                                       5603, 7781, 3568, 3668,
                                                                       6023, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 8453, 0, 3, 7781,
                                                                       5693, 7907, 3668, 3768,
                                                                       6173, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 8663, 0, 3, 8033,
                                                                       5873, 8243, 3968, 4118,
                                                                       6323, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 8978, 0, 3, 8243,
                                                                       6023, 8453, 4118, 4268,
                                                                       6548, ncols, gamma, p,
                                                                       q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 9293, 0, 3, 8663,
                                                                       6323, 8978, 4568, 4778,
                                                                       6773, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 9734, 9293, 441, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 10175, 9734, 21, 1, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 10175, 11, nmax);
    }

    for (size_t m = 0; m < 121; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
