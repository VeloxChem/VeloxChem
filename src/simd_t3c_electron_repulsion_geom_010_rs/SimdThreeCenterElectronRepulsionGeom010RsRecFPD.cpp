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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecFPD.hpp"

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
#include "SimdGeometryP1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferGeom010XDD.hpp"
#include "SimdTransferGeom010XDP.hpp"
#include "SimdTransferGeom010XFP.hpp"
#include "SimdTransferGeom010XPD.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPP.hpp"
#include "SimdTransferGeom010YDD.hpp"
#include "SimdTransferGeom010YDP.hpp"
#include "SimdTransferGeom010YFP.hpp"
#include "SimdTransferGeom010YPD.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPP.hpp"
#include "SimdTransferGeom010ZDD.hpp"
#include "SimdTransferGeom010ZDP.hpp"
#include "SimdTransferGeom010ZFP.hpp"
#include "SimdTransferGeom010ZPD.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPP.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransferPP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_fpd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_fpd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 11468, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 630 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 11468, 3976, 2587, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 7,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 15, 3, 7,
                                                             ncols, fj, i * nprim_b + j, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 484, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 487, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 490, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 493, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 496, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 499, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 502, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 505, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 508, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 511, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 514, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 517, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 520, 3, 9, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 529, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 538, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 547, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 556, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 565, 3, 18, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 574, 3, 19, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 583, 3, 20, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 592, 3, 21, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 601, 3, 22, 63,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 610, 3, 30, 78,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 628, 3, 33, 84,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 646, 3, 36, 90,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 664, 3, 39, 96,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 682, 3, 51, 114,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 700, 3, 54, 120,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 718, 3, 57, 126,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 736, 3, 60, 132,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 754, 3, 78, 158,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 784, 3, 84, 168,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 814, 3, 90, 178,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 844, 3, 114, 208,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 874, 3, 120, 218,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 904, 3, 126, 228,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 934, 3, 158, 268,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 979, 3, 168, 283,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1024, 3, 208, 328,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1069, 3, 218, 343,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1114, 3, 268, 400,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1177, 3, 328, 463,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1240, 3, 7, 8,
                                                                       484, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1246, 3, 8, 9,
                                                                       487, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1252, 3, 9, 10,
                                                                       490, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1258, 3, 10, 11,
                                                                       493, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1264, 3, 11, 12,
                                                                       496, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1270, 3, 12, 13,
                                                                       499, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1276, 3, 16, 17,
                                                                       502, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1282, 3, 17, 18,
                                                                       505, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1288, 3, 18, 19,
                                                                       508, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1294, 3, 19, 20,
                                                                       511, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1300, 3, 20, 21,
                                                                       514, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1306, 3, 21, 22,
                                                                       517, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1312, 0, 3, 1240,
                                                                       484, 1246, 520, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1330, 0, 3, 1246,
                                                                       487, 1252, 529, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1348, 0, 3, 1252,
                                                                       490, 1258, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1366, 0, 3, 1258,
                                                                       493, 1264, 547, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1384, 0, 3, 1264,
                                                                       496, 1270, 556, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1402, 0, 3, 1276,
                                                                       502, 1282, 565, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1420, 0, 3, 1282,
                                                                       505, 1288, 574, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1438, 0, 3, 1288,
                                                                       508, 1294, 583, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1456, 0, 3, 1294,
                                                                       511, 1300, 592, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1474, 0, 3, 1300,
                                                                       514, 1306, 601, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 1312,
                                                                       520, 1330, 66, 72, 610,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1528, 0, 3, 1330,
                                                                       529, 1348, 72, 78, 628,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1564, 0, 3, 1348,
                                                                       538, 1366, 78, 84, 646,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1600, 0, 3, 1366,
                                                                       547, 1384, 84, 90, 664,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1636, 0, 3, 1402,
                                                                       565, 1420, 102, 108, 682,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1672, 0, 3, 1420,
                                                                       574, 1438, 108, 114, 700,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1708, 0, 3, 1438,
                                                                       583, 1456, 114, 120, 718,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 1456,
                                                                       592, 1474, 120, 126, 736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1780, 0, 3, 1492,
                                                                       610, 1528, 138, 148, 754,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1840, 0, 3, 1528,
                                                                       628, 1564, 148, 158, 784,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1900, 0, 3, 1564,
                                                                       646, 1600, 158, 168, 814,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1960, 0, 3, 1636,
                                                                       682, 1672, 188, 198, 844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2020, 0, 3, 1672,
                                                                       700, 1708, 198, 208, 874,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2080, 0, 3, 1708,
                                                                       718, 1744, 208, 218, 904,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2140, 0, 3, 1780,
                                                                       754, 1840, 238, 253, 934,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2230, 0, 3, 1840,
                                                                       784, 1900, 253, 268, 979,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2320, 0, 3, 1960,
                                                                       844, 2020, 298, 313, 1024,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2410, 0, 3, 2020,
                                                                       874, 2080, 313, 328, 1069,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2500, 0, 3, 2140,
                                                                       934, 2230, 358, 379, 1114,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2626, 0, 3, 2320,
                                                                       1024, 2410, 421, 442,
                                                                       1177, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_p_x(buffer, 2752, 1240, 1492, 1, 6, ncols, beta);

                    simdgeo::geom_p_y(buffer, 2770, 1240, 1492, 1, 6, ncols, beta);

                    simdgeo::geom_p_z(buffer, 2788, 1240, 1492, 1, 6, ncols, beta);

                    simdgeo::geom_p_x(buffer, 2806, 1276, 1636, 1, 6, ncols, beta);

                    simdgeo::geom_p_y(buffer, 2824, 1276, 1636, 1, 6, ncols, beta);

                    simdgeo::geom_p_z(buffer, 2842, 1276, 1636, 1, 6, ncols, beta);

                    simdgeo::geom_d_x(buffer, 2860, 1312, 1780, 1, 6, ncols, beta);

                    simdgeo::geom_d_y(buffer, 2896, 1312, 1780, 1, 6, ncols, beta);

                    simdgeo::geom_d_z(buffer, 2932, 1312, 1780, 1, 6, ncols, beta);

                    simdgeo::geom_d_x(buffer, 2968, 1402, 1960, 1, 6, ncols, beta);

                    simdgeo::geom_d_y(buffer, 3004, 1402, 1960, 1, 6, ncols, beta);

                    simdgeo::geom_d_z(buffer, 3040, 1402, 1960, 1, 6, ncols, beta);

                    simdgeo::geom_f_x(buffer, 3076, 1492, 2140, 1, 6, ncols, beta);

                    simdgeo::geom_f_y(buffer, 3136, 1492, 2140, 1, 6, ncols, beta);

                    simdgeo::geom_f_z(buffer, 3196, 1492, 2140, 1, 6, ncols, beta);

                    simdgeo::geom_f_x(buffer, 3256, 1636, 2320, 1, 6, ncols, beta);

                    simdgeo::geom_f_y(buffer, 3316, 1636, 2320, 1, 6, ncols, beta);

                    simdgeo::geom_f_z(buffer, 3376, 1636, 2320, 1, 6, ncols, beta);

                    simdgeo::geom_g_x(buffer, 3436, 1780, 2500, 1, 6, ncols, beta);

                    simdgeo::geom_g_y(buffer, 3526, 1780, 2500, 1, 6, ncols, beta);

                    simdgeo::geom_g_z(buffer, 3616, 1780, 2500, 1, 6, ncols, beta);

                    simdgeo::geom_g_x(buffer, 3706, 1960, 2626, 1, 6, ncols, beta);

                    simdgeo::geom_g_y(buffer, 3796, 1960, 2626, 1, 6, ncols, beta);

                    simdgeo::geom_g_z(buffer, 3886, 1960, 2626, 1, 6, ncols, beta);

                    simdfunc::contract_primitives(buffer, 3976, 2752, 18, ncols);

                    simdfunc::contract_primitives(buffer, 4009, 2770, 18, ncols);

                    simdfunc::contract_primitives(buffer, 4042, 2788, 18, ncols);

                    simdfunc::contract_primitives(buffer, 4075, 1312, 18, ncols);

                    simdfunc::contract_primitives(buffer, 4108, 2806, 18, ncols);

                    simdfunc::contract_primitives(buffer, 4141, 2824, 18, ncols);

                    simdfunc::contract_primitives(buffer, 4174, 2842, 18, ncols);

                    simdfunc::contract_primitives(buffer, 4207, 1402, 18, ncols);

                    simdfunc::contract_primitives(buffer, 4240, 2860, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4306, 2896, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4372, 2932, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4438, 1492, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4504, 2968, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4570, 3004, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4636, 3040, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4702, 1636, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4768, 3076, 60, ncols);

                    simdfunc::contract_primitives(buffer, 4878, 3136, 60, ncols);

                    simdfunc::contract_primitives(buffer, 4988, 3196, 60, ncols);

                    simdfunc::contract_primitives(buffer, 5098, 1780, 60, ncols);

                    simdfunc::contract_primitives(buffer, 5208, 3256, 60, ncols);

                    simdfunc::contract_primitives(buffer, 5318, 3316, 60, ncols);

                    simdfunc::contract_primitives(buffer, 5428, 3376, 60, ncols);

                    simdfunc::contract_primitives(buffer, 5538, 1960, 60, ncols);

                    simdfunc::contract_primitives(buffer, 5648, 3436, 90, ncols);

                    simdfunc::contract_primitives(buffer, 5813, 3526, 90, ncols);

                    simdfunc::contract_primitives(buffer, 5978, 3616, 90, ncols);

                    simdfunc::contract_primitives(buffer, 6143, 3706, 90, ncols);

                    simdfunc::contract_primitives(buffer, 6308, 3796, 90, ncols);

                    simdfunc::contract_primitives(buffer, 6473, 3886, 90, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 3994, 3976, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4027, 4009, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4060, 4042, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4093, 4075, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4126, 4108, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4159, 4141, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4192, 4174, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4225, 4207, 3, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4276, 4240, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4342, 4306, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4408, 4372, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4474, 4438, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4540, 4504, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4606, 4570, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4672, 4636, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4738, 4702, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4828, 4768, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4938, 4878, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5048, 4988, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5158, 5098, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5268, 5208, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5378, 5318, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5488, 5428, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5598, 5538, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5738, 5648, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5903, 5813, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 6068, 5978, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 6233, 6143, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 6398, 6308, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 6563, 6473, 15, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 6638, 3994, 4093, 4276, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 6683, 4027, 4093, 4342, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 6728, 4060, 4093, 4408, 5, nmax);

        simdtrf::compute_hrr_pp(buffer, coordinates, 6773, 4093, 4474, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 6818, 4126, 4225, 4540, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 6863, 4159, 4225, 4606, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 6908, 4192, 4225, 4672, 5, nmax);

        simdtrf::compute_hrr_pp(buffer, coordinates, 6953, 4225, 4738, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 6998, 4276, 4474, 4828, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 7088, 4342, 4474, 4938, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 7178, 4408, 4474, 5048, 5, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 7268, 4474, 5158, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 7358, 4540, 4738, 5268, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 7448, 4606, 4738, 5378, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 7538, 4672, 4738, 5488, 5, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 7628, 4738, 5598, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 7718, 4828, 5158, 5738, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 7868, 4938, 5158, 5903, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 8018, 5048, 5158, 6068, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 8168, 5268, 5598, 6233, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 8318, 5378, 5598, 6398, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 8468, 5488, 5598, 6563, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dp_out_of_second(buffer, coordinates, 8618, 6638, 6773,
                                                        6998, 5, nmax);

        simdtrf::compute_hrr_geom_010y_dp_out_of_second(buffer, coordinates, 8708, 6683, 6773,
                                                        7088, 5, nmax);

        simdtrf::compute_hrr_geom_010z_dp_out_of_second(buffer, coordinates, 8798, 6728, 6773,
                                                        7178, 5, nmax);

        simdtrf::compute_hrr_dp_out_of_second(buffer, coordinates, 8888, 6773, 7268, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dp_out_of_second(buffer, coordinates, 8978, 6818, 6953,
                                                        7358, 5, nmax);

        simdtrf::compute_hrr_geom_010y_dp_out_of_second(buffer, coordinates, 9068, 6863, 6953,
                                                        7448, 5, nmax);

        simdtrf::compute_hrr_geom_010z_dp_out_of_second(buffer, coordinates, 9158, 6908, 6953,
                                                        7538, 5, nmax);

        simdtrf::compute_hrr_dp_out_of_second(buffer, coordinates, 9248, 6953, 7628, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 9338, 6998, 7268, 7718, 5, nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 9518, 7088, 7268, 7868, 5, nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 9698, 7178, 7268, 8018, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 9878, 7358, 7628, 8168, 5, nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 10058, 7448, 7628, 8318, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 10238, 7538, 7628, 8468, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fp_out_of_second(buffer, coordinates, 10418, 8618, 8888,
                                                        9338, 5, nmax);

        simdtrf::compute_hrr_geom_010y_fp_out_of_second(buffer, coordinates, 10568, 8708, 8888,
                                                        9518, 5, nmax);

        simdtrf::compute_hrr_geom_010z_fp_out_of_second(buffer, coordinates, 10718, 8798, 8888,
                                                        9698, 5, nmax);

        simdtrf::compute_hrr_geom_010x_fp_out_of_second(buffer, coordinates, 10868, 8978, 9248,
                                                        9878, 5, nmax);

        simdtrf::compute_hrr_geom_010y_fp_out_of_second(buffer, coordinates, 11018, 9068, 9248,
                                                        10058, 5, nmax);

        simdtrf::compute_hrr_geom_010z_fp_out_of_second(buffer, coordinates, 11168, 9158, 9248,
                                                        10238, 5, nmax);

        simdtrf::transform_p_inner(buffer, 11318, 10868, 10, 5, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 11318, 15, nmax);

        simdtrf::transform_p_inner(buffer, 11318, 11018, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 105 * nvalues + n * npairs, nvalues, buffer, 11318,
                                   15, nmax);

        simdtrf::transform_p_inner(buffer, 11318, 11168, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 210 * nvalues + n * npairs, nvalues, buffer, 11318,
                                   15, nmax);

        simdtrf::transform_p_inner(buffer, 11318, 10418, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 11318,
                                   15, nmax);

        simdtrf::transform_p_inner(buffer, 11318, 10568, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 420 * nvalues + n * npairs, nvalues, buffer, 11318,
                                   15, nmax);

        simdtrf::transform_p_inner(buffer, 11318, 10718, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 525 * nvalues + n * npairs, nvalues, buffer, 11318,
                                   15, nmax);
    }

    for (size_t m = 0; m < 630; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
