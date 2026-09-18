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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecDDF.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdGeometryD1.hpp"
#include "SimdGeometryF1.hpp"
#include "SimdGeometryG1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferGeom010XDD.hpp"
#include "SimdTransferGeom010XPD.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010YDD.hpp"
#include "SimdTransferGeom010YPD.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010ZDD.hpp"
#include "SimdTransferGeom010ZPD.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_ddf_three_center_electron_repulsion(double               *values,
                                                        const size_t          npairs,
                                                        const size_t          natoms,
                                                        const CBasisFunction &a_function,
                                                        const CBasisFunction &b_function,
                                                        const CBasisFunction &c_function,
                                                        const CSimdMatrix    &coordinates,
                                                        const CSimdMatrix    &c_coordinates,
                                                        CSimdMatrix          &buffer,
                                                        const double          omega,
                                                        const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_geom_010_ddf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    // NOTE: a derivative screens with the integral's own bound, on purpose. It
    // is not a bound on the derivative -- that is larger by roughly 2 alpha R,
    // the relation reaching one shell higher -- and tightening it here would be
    // the wrong repair. A screened Fock build defines an energy in which the
    // dropped pairs contribute exactly zero, and the derivative of that energy
    // is the derivative screened the same way; a tighter bound would add forces
    // from pairs the energy never counted. One threshold controls both errors,
    // so tightening it in the Fock build tightens the gradient with it.

    const auto dimensions = simdfunc::make_column_dimensions(
        a_function, b_function, c_function, npairs, coordinates,
        screenfunc::three_center_electron_repulsion_primitive_bound,
        threshold / static_cast<double>(nprims));

    const auto nmax = simdfunc::prepare_buffer(buffer, 15464, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1050 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    simdfunc::compute_pair_exponents(a_function, b_function, coordinates, nmax);

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 15464, 7768, 3601, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto beta = b_exps[j];

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pb(buffer, coordinates, 0, nmax, fb);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 3, n, nmax, fc);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_t3c_erf_boys_function(buffer, coordinates, 6, 3, {1, 2, 3,
                                                            4, 5, 6, 7, 8}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 15, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 60, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 63, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 66, 0, 3, 7, 8,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 72, 0, 3, 8, 9,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 78, 0, 3, 9, 10,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 84, 0, 3, 10, 11,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 90, 0, 3, 11, 12,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 96, 0, 3, 12, 13,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 102, 0, 3, 16, 17,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 108, 0, 3, 17, 18,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 114, 0, 3, 18, 19,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 120, 0, 3, 19, 20,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 126, 0, 3, 20, 21,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 132, 0, 3, 21, 22,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 138, 0, 3, 24, 27,
                                                                       66, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 148, 0, 3, 27, 30,
                                                                       72, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 158, 0, 3, 30, 33,
                                                                       78, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 168, 0, 3, 33, 36,
                                                                       84, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 178, 0, 3, 36, 39,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 188, 0, 3, 45, 48,
                                                                       102, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 198, 0, 3, 48, 51,
                                                                       108, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 208, 0, 3, 51, 54,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 218, 0, 3, 54, 57,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 57, 60,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 66, 72,
                                                                       138, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 253, 0, 3, 72, 78,
                                                                       148, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 78, 84,
                                                                       158, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 283, 0, 3, 84, 90,
                                                                       168, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 102,
                                                                       108, 188, 198, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 313, 0, 3, 108,
                                                                       114, 198, 208, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 114,
                                                                       120, 208, 218, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 343, 0, 3, 120,
                                                                       126, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 138,
                                                                       148, 238, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 379, 0, 3, 148,
                                                                       158, 253, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 400, 0, 3, 158,
                                                                       168, 268, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 421, 0, 3, 188,
                                                                       198, 298, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 442, 0, 3, 198,
                                                                       208, 313, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 463, 0, 3, 208,
                                                                       218, 328, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 484, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 487, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 490, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 493, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 496, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 499, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 502, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 505, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 508, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 511, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 514, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 517, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 520, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 523, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 526, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 529, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 532, 3, 9, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 541, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 550, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 559, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 568, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 577, 3, 18, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 586, 3, 19, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 595, 3, 20, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 604, 3, 21, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 613, 3, 22, 63,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 622, 3, 24, 66,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 640, 3, 27, 72,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 658, 3, 30, 78,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 676, 3, 33, 84,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 694, 3, 36, 90,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 712, 3, 39, 96,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 730, 3, 45, 102,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 748, 3, 48, 108,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 766, 3, 51, 114,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 784, 3, 54, 120,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 802, 3, 57, 126,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 820, 3, 60, 132,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 838, 3, 66, 138,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 868, 3, 72, 148,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 898, 3, 78, 158,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 928, 3, 84, 168,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 958, 3, 90, 178,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 988, 3, 102, 188,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1018, 3, 108, 198,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1048, 3, 114, 208,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1078, 3, 120, 218,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1108, 3, 126, 228,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1138, 3, 138, 238,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1183, 3, 148, 253,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1228, 3, 158, 268,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1273, 3, 168, 283,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1318, 3, 188, 298,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1363, 3, 198, 313,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1408, 3, 208, 328,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1453, 3, 218, 343,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1498, 3, 238, 358,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1561, 3, 253, 379,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1624, 3, 268, 400,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1687, 3, 298, 421,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1750, 3, 313, 442,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1813, 3, 328, 463,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1876, 3, 7, 8,
                                                                       490, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1882, 3, 8, 9,
                                                                       493, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1888, 3, 9, 10,
                                                                       496, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1894, 3, 10, 11,
                                                                       499, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1900, 3, 11, 12,
                                                                       502, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1906, 3, 12, 13,
                                                                       505, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1912, 3, 16, 17,
                                                                       514, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1918, 3, 17, 18,
                                                                       517, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1924, 3, 18, 19,
                                                                       520, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1930, 3, 19, 20,
                                                                       523, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1936, 3, 20, 21,
                                                                       526, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1942, 3, 21, 22,
                                                                       529, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1948, 0, 3, 1876,
                                                                       490, 1882, 532, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1966, 0, 3, 1882,
                                                                       493, 1888, 541, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1984, 0, 3, 1888,
                                                                       496, 1894, 550, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2002, 0, 3, 1894,
                                                                       499, 1900, 559, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2020, 0, 3, 1900,
                                                                       502, 1906, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2038, 0, 3, 1912,
                                                                       514, 1918, 577, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2056, 0, 3, 1918,
                                                                       517, 1924, 586, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2074, 0, 3, 1924,
                                                                       520, 1930, 595, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2092, 0, 3, 1930,
                                                                       523, 1936, 604, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2110, 0, 3, 1936,
                                                                       526, 1942, 613, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2128, 0, 3, 1948,
                                                                       532, 1966, 66, 72, 658,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2164, 0, 3, 1966,
                                                                       541, 1984, 72, 78, 676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2200, 0, 3, 1984,
                                                                       550, 2002, 78, 84, 694,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2236, 0, 3, 2002,
                                                                       559, 2020, 84, 90, 712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2272, 0, 3, 2038,
                                                                       577, 2056, 102, 108, 766,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2308, 0, 3, 2056,
                                                                       586, 2074, 108, 114, 784,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2344, 0, 3, 2074,
                                                                       595, 2092, 114, 120, 802,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2380, 0, 3, 2092,
                                                                       604, 2110, 120, 126, 820,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2416, 0, 3, 2128,
                                                                       658, 2164, 138, 148, 898,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2476, 0, 3, 2164,
                                                                       676, 2200, 148, 158, 928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2536, 0, 3, 2200,
                                                                       694, 2236, 158, 168, 958,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2596, 0, 3, 2272,
                                                                       766, 2308, 188, 198, 1048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2656, 0, 3, 2308,
                                                                       784, 2344, 198, 208, 1078,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2716, 0, 3, 2344,
                                                                       802, 2380, 208, 218, 1108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2776, 0, 3, 2416,
                                                                       898, 2476, 238, 253, 1228,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2866, 0, 3, 2476,
                                                                       928, 2536, 253, 268, 1273,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2956, 0, 3, 2596,
                                                                       1048, 2656, 298, 313,
                                                                       1408, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3046, 0, 3, 2656,
                                                                       1078, 2716, 313, 328,
                                                                       1453, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3136, 0, 3, 2776,
                                                                       1228, 2866, 358, 379,
                                                                       1624, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3262, 0, 3, 2956,
                                                                       1408, 3046, 421, 442,
                                                                       1813, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3388, 3, 484, 487,
                                                                       1876, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3398, 3, 487, 490,
                                                                       1882, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3408, 3, 490, 493,
                                                                       1888, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3418, 3, 493, 496,
                                                                       1894, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3428, 3, 496, 499,
                                                                       1900, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3438, 3, 499, 502,
                                                                       1906, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3448, 3, 508, 511,
                                                                       1912, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3458, 3, 511, 514,
                                                                       1918, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3468, 3, 514, 517,
                                                                       1924, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3478, 3, 517, 520,
                                                                       1930, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3488, 3, 520, 523,
                                                                       1936, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3498, 3, 523, 526,
                                                                       1942, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3508, 0, 3, 3388,
                                                                       1876, 3398, 1948, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3538, 0, 3, 3398,
                                                                       1882, 3408, 1966, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3568, 0, 3, 3408,
                                                                       1888, 3418, 1984, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3598, 0, 3, 3418,
                                                                       1894, 3428, 2002, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3628, 0, 3, 3428,
                                                                       1900, 3438, 2020, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3658, 0, 3, 3448,
                                                                       1912, 3458, 2038, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3688, 0, 3, 3458,
                                                                       1918, 3468, 2056, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3718, 0, 3, 3468,
                                                                       1924, 3478, 2074, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3748, 0, 3, 3478,
                                                                       1930, 3488, 2092, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3778, 0, 3, 3488,
                                                                       1936, 3498, 2110, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3808, 0, 3, 3508,
                                                                       1948, 3538, 622, 640,
                                                                       2128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3868, 0, 3, 3538,
                                                                       1966, 3568, 640, 658,
                                                                       2164, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3928, 0, 3, 3568,
                                                                       1984, 3598, 658, 676,
                                                                       2200, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3988, 0, 3, 3598,
                                                                       2002, 3628, 676, 694,
                                                                       2236, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4048, 0, 3, 3658,
                                                                       2038, 3688, 730, 748,
                                                                       2272, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4108, 0, 3, 3688,
                                                                       2056, 3718, 748, 766,
                                                                       2308, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4168, 0, 3, 3718,
                                                                       2074, 3748, 766, 784,
                                                                       2344, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4228, 0, 3, 3748,
                                                                       2092, 3778, 784, 802,
                                                                       2380, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4288, 0, 3, 3808,
                                                                       2128, 3868, 838, 868,
                                                                       2416, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4388, 0, 3, 3868,
                                                                       2164, 3928, 868, 898,
                                                                       2476, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4488, 0, 3, 3928,
                                                                       2200, 3988, 898, 928,
                                                                       2536, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4588, 0, 3, 4048,
                                                                       2272, 4108, 988, 1018,
                                                                       2596, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4688, 0, 3, 4108,
                                                                       2308, 4168, 1018, 1048,
                                                                       2656, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4788, 0, 3, 4168,
                                                                       2344, 4228, 1048, 1078,
                                                                       2716, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 4888, 0, 3, 4288,
                                                                       2416, 4388, 1138, 1183,
                                                                       2776, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5038, 0, 3, 4388,
                                                                       2476, 4488, 1183, 1228,
                                                                       2866, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5188, 0, 3, 4588,
                                                                       2596, 4688, 1318, 1363,
                                                                       2956, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5338, 0, 3, 4688,
                                                                       2656, 4788, 1363, 1408,
                                                                       3046, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 5488, 0, 3, 4888,
                                                                       2776, 5038, 1498, 1561,
                                                                       3136, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 5698, 0, 3, 5188,
                                                                       2956, 5338, 1687, 1750,
                                                                       3262, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_d_x(buffer, 5908, 3508, 4288, 1, 10, ncols, beta);

                    simdgeo::geom_d_y(buffer, 5968, 3508, 4288, 1, 10, ncols, beta);

                    simdgeo::geom_d_z(buffer, 6028, 3508, 4288, 1, 10, ncols, beta);

                    simdgeo::geom_d_x(buffer, 6088, 3658, 4588, 1, 10, ncols, beta);

                    simdgeo::geom_d_y(buffer, 6148, 3658, 4588, 1, 10, ncols, beta);

                    simdgeo::geom_d_z(buffer, 6208, 3658, 4588, 1, 10, ncols, beta);

                    simdgeo::geom_f_x(buffer, 6268, 3808, 4888, 1, 10, ncols, beta);

                    simdgeo::geom_f_y(buffer, 6368, 3808, 4888, 1, 10, ncols, beta);

                    simdgeo::geom_f_z(buffer, 6468, 3808, 4888, 1, 10, ncols, beta);

                    simdgeo::geom_f_x(buffer, 6568, 4048, 5188, 1, 10, ncols, beta);

                    simdgeo::geom_f_y(buffer, 6668, 4048, 5188, 1, 10, ncols, beta);

                    simdgeo::geom_f_z(buffer, 6768, 4048, 5188, 1, 10, ncols, beta);

                    simdgeo::geom_g_x(buffer, 6868, 4288, 5488, 1, 10, ncols, beta);

                    simdgeo::geom_g_y(buffer, 7018, 4288, 5488, 1, 10, ncols, beta);

                    simdgeo::geom_g_z(buffer, 7168, 4288, 5488, 1, 10, ncols, beta);

                    simdgeo::geom_g_x(buffer, 7318, 4588, 5698, 1, 10, ncols, beta);

                    simdgeo::geom_g_y(buffer, 7468, 4588, 5698, 1, 10, ncols, beta);

                    simdgeo::geom_g_z(buffer, 7618, 4588, 5698, 1, 10, ncols, beta);

                    simdfunc::contract_primitives(buffer, 7768, 5908, 60, ncols);

                    simdfunc::contract_primitives(buffer, 7870, 5968, 60, ncols);

                    simdfunc::contract_primitives(buffer, 7972, 6028, 60, ncols);

                    simdfunc::contract_primitives(buffer, 8074, 3808, 60, ncols);

                    simdfunc::contract_primitives(buffer, 8176, 6088, 60, ncols);

                    simdfunc::contract_primitives(buffer, 8278, 6148, 60, ncols);

                    simdfunc::contract_primitives(buffer, 8380, 6208, 60, ncols);

                    simdfunc::contract_primitives(buffer, 8482, 4048, 60, ncols);

                    simdfunc::contract_primitives(buffer, 8584, 6268, 100, ncols);

                    simdfunc::contract_primitives(buffer, 8754, 6368, 100, ncols);

                    simdfunc::contract_primitives(buffer, 8924, 6468, 100, ncols);

                    simdfunc::contract_primitives(buffer, 9094, 4288, 100, ncols);

                    simdfunc::contract_primitives(buffer, 9264, 6568, 100, ncols);

                    simdfunc::contract_primitives(buffer, 9434, 6668, 100, ncols);

                    simdfunc::contract_primitives(buffer, 9604, 6768, 100, ncols);

                    simdfunc::contract_primitives(buffer, 9774, 4588, 100, ncols);

                    simdfunc::contract_primitives(buffer, 9944, 6868, 150, ncols);

                    simdfunc::contract_primitives(buffer, 10199, 7018, 150, ncols);

                    simdfunc::contract_primitives(buffer, 10454, 7168, 150, ncols);

                    simdfunc::contract_primitives(buffer, 10709, 7318, 150, ncols);

                    simdfunc::contract_primitives(buffer, 10964, 7468, 150, ncols);

                    simdfunc::contract_primitives(buffer, 11219, 7618, 150, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 7828, 7768, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 7930, 7870, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8032, 7972, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8134, 8074, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8236, 8176, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8338, 8278, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8440, 8380, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8542, 8482, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8684, 8584, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8854, 8754, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 9024, 8924, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 9194, 9094, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 9364, 9264, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 9534, 9434, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 9704, 9604, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 9874, 9774, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10094, 9944, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10349, 10199, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10604, 10454, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10859, 10709, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 11114, 10964, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 11369, 11219, 15, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 11474, 7828, 8134, 8684, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 11600, 7930, 8134, 8854, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 11726, 8032, 8134, 9024, 7,
                                          nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 11852, 8134, 9194, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 11978, 8236, 8542, 9364, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 12104, 8338, 8542, 9534, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 12230, 8440, 8542, 9704, 7,
                                          nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 12356, 8542, 9874, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 12482, 8684, 9194, 10094, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 12692, 8854, 9194, 10349, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 12902, 9024, 9194, 10604, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 13112, 9364, 9874, 10859, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 13322, 9534, 9874, 11114, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 13532, 9704, 9874, 11369, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 13742, 11474, 11852, 12482, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 13994, 11600, 11852, 12692, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 14246, 11726, 11852, 12902, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 14498, 11978, 12356, 13112, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 14750, 12104, 12356, 13322, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 15002, 12230, 12356, 13532, 7,
                                          nmax);

        simdtrf::transform_d_inner(buffer, 15254, 14498, 6, 7, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 15254, 35, nmax);

        simdtrf::transform_d_inner(buffer, 15254, 14750, 6, 7, nmax);

        simdtrf::transform_d_outer(values + 175 * nvalues + n * npairs, nvalues, buffer, 15254,
                                   35, nmax);

        simdtrf::transform_d_inner(buffer, 15254, 15002, 6, 7, nmax);

        simdtrf::transform_d_outer(values + 350 * nvalues + n * npairs, nvalues, buffer, 15254,
                                   35, nmax);

        simdtrf::transform_d_inner(buffer, 15254, 13742, 6, 7, nmax);

        simdtrf::transform_d_outer(values + 525 * nvalues + n * npairs, nvalues, buffer, 15254,
                                   35, nmax);

        simdtrf::transform_d_inner(buffer, 15254, 13994, 6, 7, nmax);

        simdtrf::transform_d_outer(values + 700 * nvalues + n * npairs, nvalues, buffer, 15254,
                                   35, nmax);

        simdtrf::transform_d_inner(buffer, 15254, 14246, 6, 7, nmax);

        simdtrf::transform_d_outer(values + 875 * nvalues + n * npairs, nvalues, buffer, 15254,
                                   35, nmax);
    }

    for (size_t m = 0; m < 1050; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
