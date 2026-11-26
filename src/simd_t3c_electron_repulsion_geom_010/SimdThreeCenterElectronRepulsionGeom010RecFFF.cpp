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


#include "SimdThreeCenterElectronRepulsionGeom010RecFFF.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdGeometryF1.hpp"
#include "SimdGeometryG1.hpp"
#include "SimdGeometryH1.hpp"
#include "SimdGeometryI1.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDF.hpp"
#include "SimdTransferGeom010XDF.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XFF.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010YDF.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YFF.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010ZDF.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZFF.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransformF.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_fff_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_fff_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 23908, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1029 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 23908, 9769, 4360, dimensions);

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 6, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 17, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 20, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 44, 0, 3, 7, 8,
                                                                       17, 20, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 50, 0, 3, 8, 9,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 56, 0, 3, 9, 10,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 62, 0, 3, 10, 11,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 68, 0, 3, 11, 12,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 74, 0, 3, 12, 13,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 80, 0, 3, 13, 14,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 86, 0, 3, 14, 15,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 92, 0, 3, 17, 20,
                                                                       44, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 102, 0, 3, 20, 23,
                                                                       50, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 112, 0, 3, 23, 26,
                                                                       56, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 122, 0, 3, 26, 29,
                                                                       62, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 132, 0, 3, 29, 32,
                                                                       68, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 142, 0, 3, 32, 35,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 152, 0, 3, 35, 38,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 162, 0, 3, 44, 50,
                                                                       92, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 177, 0, 3, 50, 56,
                                                                       102, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 192, 0, 3, 56, 62,
                                                                       112, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 207, 0, 3, 62, 68,
                                                                       122, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 222, 0, 3, 68, 74,
                                                                       132, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 237, 0, 3, 74, 80,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 252, 0, 3, 92,
                                                                       102, 162, 177, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 102,
                                                                       112, 177, 192, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 294, 0, 3, 112,
                                                                       122, 192, 207, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 315, 0, 3, 122,
                                                                       132, 207, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 336, 0, 3, 132,
                                                                       142, 222, 237, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 357, 0, 3, 162,
                                                                       177, 252, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 385, 0, 3, 177,
                                                                       192, 273, 294, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 413, 0, 3, 192,
                                                                       207, 294, 315, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 441, 0, 3, 207,
                                                                       222, 315, 336, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 469, 0, 3, 252,
                                                                       273, 357, 385, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 505, 0, 3, 273,
                                                                       294, 385, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 541, 0, 3, 294,
                                                                       315, 413, 441, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 577, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 580, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 583, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 586, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 589, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 592, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 595, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 598, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 601, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 604, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 607, 3, 9, 23,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 616, 3, 10, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 625, 3, 11, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 634, 3, 12, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 643, 3, 13, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 652, 3, 14, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 661, 3, 15, 41,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 670, 3, 17, 44,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 688, 3, 20, 50,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 706, 3, 23, 56,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 724, 3, 26, 62,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 742, 3, 29, 68,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 760, 3, 32, 74,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 778, 3, 35, 80,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 796, 3, 38, 86,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 814, 3, 44, 92,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 844, 3, 50, 102,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 874, 3, 56, 112,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 904, 3, 62, 122,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 934, 3, 68, 132,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 964, 3, 74, 142,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 994, 3, 80, 152,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1024, 3, 92, 162,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1069, 3, 102, 177,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1114, 3, 112, 192,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1159, 3, 122, 207,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1204, 3, 132, 222,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1249, 3, 142, 237,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1294, 3, 162, 252,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1357, 3, 177, 273,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1420, 3, 192, 294,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1483, 3, 207, 315,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1546, 3, 222, 336,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1609, 3, 252, 357,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1693, 3, 273, 385,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1777, 3, 294, 413,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1861, 3, 315, 441,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 1945, 3, 357, 469,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 2053, 3, 385, 505,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 2161, 3, 413, 541,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2269, 3, 7, 8,
                                                                       583, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2275, 3, 8, 9,
                                                                       586, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2281, 3, 9, 10,
                                                                       589, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2287, 3, 10, 11,
                                                                       592, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2293, 3, 11, 12,
                                                                       595, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2299, 3, 12, 13,
                                                                       598, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2305, 3, 13, 14,
                                                                       601, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2311, 3, 14, 15,
                                                                       604, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2317, 0, 3, 2269,
                                                                       583, 2275, 607, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2335, 0, 3, 2275,
                                                                       586, 2281, 616, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2353, 0, 3, 2281,
                                                                       589, 2287, 625, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2371, 0, 3, 2287,
                                                                       592, 2293, 634, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2389, 0, 3, 2293,
                                                                       595, 2299, 643, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2407, 0, 3, 2299,
                                                                       598, 2305, 652, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2425, 0, 3, 2305,
                                                                       601, 2311, 661, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2443, 0, 3, 2317,
                                                                       607, 2335, 44, 50, 706,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2479, 0, 3, 2335,
                                                                       616, 2353, 50, 56, 724,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2515, 0, 3, 2353,
                                                                       625, 2371, 56, 62, 742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2551, 0, 3, 2371,
                                                                       634, 2389, 62, 68, 760,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2587, 0, 3, 2389,
                                                                       643, 2407, 68, 74, 778,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2623, 0, 3, 2407,
                                                                       652, 2425, 74, 80, 796,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2659, 0, 3, 2443,
                                                                       706, 2479, 92, 102, 874,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2719, 0, 3, 2479,
                                                                       724, 2515, 102, 112, 904,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2779, 0, 3, 2515,
                                                                       742, 2551, 112, 122, 934,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2839, 0, 3, 2551,
                                                                       760, 2587, 122, 132, 964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2899, 0, 3, 2587,
                                                                       778, 2623, 132, 142, 994,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2959, 0, 3, 2659,
                                                                       874, 2719, 162, 177, 1114,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3049, 0, 3, 2719,
                                                                       904, 2779, 177, 192, 1159,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3139, 0, 3, 2779,
                                                                       934, 2839, 192, 207, 1204,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3229, 0, 3, 2839,
                                                                       964, 2899, 207, 222, 1249,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3319, 0, 3, 2959,
                                                                       1114, 3049, 252, 273,
                                                                       1420, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3445, 0, 3, 3049,
                                                                       1159, 3139, 273, 294,
                                                                       1483, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3571, 0, 3, 3139,
                                                                       1204, 3229, 294, 315,
                                                                       1546, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 3697, 0, 3, 3319,
                                                                       1420, 3445, 357, 385,
                                                                       1777, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 3865, 0, 3, 3445,
                                                                       1483, 3571, 385, 413,
                                                                       1861, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 4033, 0, 3, 3697,
                                                                       1777, 3865, 469, 505,
                                                                       2161, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4249, 3, 577, 580,
                                                                       2269, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4259, 3, 580, 583,
                                                                       2275, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4269, 3, 583, 586,
                                                                       2281, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4279, 3, 586, 589,
                                                                       2287, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4289, 3, 589, 592,
                                                                       2293, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4299, 3, 592, 595,
                                                                       2299, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4309, 3, 595, 598,
                                                                       2305, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4319, 3, 598, 601,
                                                                       2311, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4329, 0, 3, 4249,
                                                                       2269, 4259, 2317, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4359, 0, 3, 4259,
                                                                       2275, 4269, 2335, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4389, 0, 3, 4269,
                                                                       2281, 4279, 2353, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4419, 0, 3, 4279,
                                                                       2287, 4289, 2371, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4449, 0, 3, 4289,
                                                                       2293, 4299, 2389, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4479, 0, 3, 4299,
                                                                       2299, 4309, 2407, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4509, 0, 3, 4309,
                                                                       2305, 4319, 2425, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4539, 0, 3, 4329,
                                                                       2317, 4359, 670, 688,
                                                                       2443, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4599, 0, 3, 4359,
                                                                       2335, 4389, 688, 706,
                                                                       2479, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4659, 0, 3, 4389,
                                                                       2353, 4419, 706, 724,
                                                                       2515, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4719, 0, 3, 4419,
                                                                       2371, 4449, 724, 742,
                                                                       2551, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4779, 0, 3, 4449,
                                                                       2389, 4479, 742, 760,
                                                                       2587, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4839, 0, 3, 4479,
                                                                       2407, 4509, 760, 778,
                                                                       2623, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4899, 0, 3, 4539,
                                                                       2443, 4599, 814, 844,
                                                                       2659, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4999, 0, 3, 4599,
                                                                       2479, 4659, 844, 874,
                                                                       2719, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5099, 0, 3, 4659,
                                                                       2515, 4719, 874, 904,
                                                                       2779, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5199, 0, 3, 4719,
                                                                       2551, 4779, 904, 934,
                                                                       2839, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 5299, 0, 3, 4779,
                                                                       2587, 4839, 934, 964,
                                                                       2899, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5399, 0, 3, 4899,
                                                                       2659, 4999, 1024, 1069,
                                                                       2959, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5549, 0, 3, 4999,
                                                                       2719, 5099, 1069, 1114,
                                                                       3049, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5699, 0, 3, 5099,
                                                                       2779, 5199, 1114, 1159,
                                                                       3139, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 5849, 0, 3, 5199,
                                                                       2839, 5299, 1159, 1204,
                                                                       3229, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 5999, 0, 3, 5399,
                                                                       2959, 5549, 1294, 1357,
                                                                       3319, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 6209, 0, 3, 5549,
                                                                       3049, 5699, 1357, 1420,
                                                                       3445, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 6419, 0, 3, 5699,
                                                                       3139, 5849, 1420, 1483,
                                                                       3571, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 6629, 0, 3, 5999,
                                                                       3319, 6209, 1609, 1693,
                                                                       3697, ncols, gamma, p,
                                                                       q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 6909, 0, 3, 6209,
                                                                       3445, 6419, 1693, 1777,
                                                                       3865, ncols, gamma, p,
                                                                       q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 7189, 0, 3, 6629,
                                                                       3697, 6909, 1945, 2053,
                                                                       4033, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_f_x(buffer, 7549, 4539, 5399, 1, 10, ncols, beta);

                    simdgeo::geom_f_y(buffer, 7649, 4539, 5399, 1, 10, ncols, beta);

                    simdgeo::geom_f_z(buffer, 7749, 4539, 5399, 1, 10, ncols, beta);

                    simdgeo::geom_g_x(buffer, 7849, 4899, 5999, 1, 10, ncols, beta);

                    simdgeo::geom_g_y(buffer, 7999, 4899, 5999, 1, 10, ncols, beta);

                    simdgeo::geom_g_z(buffer, 8149, 4899, 5999, 1, 10, ncols, beta);

                    simdgeo::geom_h_x(buffer, 8299, 5399, 6629, 1, 10, ncols, beta);

                    simdgeo::geom_h_y(buffer, 8509, 5399, 6629, 1, 10, ncols, beta);

                    simdgeo::geom_h_z(buffer, 8719, 5399, 6629, 1, 10, ncols, beta);

                    simdgeo::geom_i_x(buffer, 8929, 5999, 7189, 1, 10, ncols, beta);

                    simdgeo::geom_i_y(buffer, 9209, 5999, 7189, 1, 10, ncols, beta);

                    simdgeo::geom_i_z(buffer, 9489, 5999, 7189, 1, 10, ncols, beta);

                    simdfunc::contract_primitives(buffer, 9769, 7549, 100, ncols);

                    simdfunc::contract_primitives(buffer, 9939, 7649, 100, ncols);

                    simdfunc::contract_primitives(buffer, 10109, 7749, 100, ncols);

                    simdfunc::contract_primitives(buffer, 10279, 4899, 100, ncols);

                    simdfunc::contract_primitives(buffer, 10449, 7849, 150, ncols);

                    simdfunc::contract_primitives(buffer, 10704, 7999, 150, ncols);

                    simdfunc::contract_primitives(buffer, 10959, 8149, 150, ncols);

                    simdfunc::contract_primitives(buffer, 11214, 5399, 150, ncols);

                    simdfunc::contract_primitives(buffer, 11469, 8299, 210, ncols);

                    simdfunc::contract_primitives(buffer, 11826, 8509, 210, ncols);

                    simdfunc::contract_primitives(buffer, 12183, 8719, 210, ncols);

                    simdfunc::contract_primitives(buffer, 12540, 5999, 210, ncols);

                    simdfunc::contract_primitives(buffer, 12897, 8929, 280, ncols);

                    simdfunc::contract_primitives(buffer, 13373, 9209, 280, ncols);

                    simdfunc::contract_primitives(buffer, 13849, 9489, 280, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 9869, 9769, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10039, 9939, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10209, 10109, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10379, 10279, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10599, 10449, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10854, 10704, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 11109, 10959, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 11364, 11214, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 11679, 11469, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 12036, 11826, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 12393, 12183, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 12750, 12540, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13177, 12897, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 13653, 13373, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 14129, 13849, 28, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 14325, 9869, 10379, 10599, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 14535, 10039, 10379, 10854, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 14745, 10209, 10379, 11109, 7,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 14955, 10379, 11364, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 15165, 10599, 11364, 11679, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 15480, 10854, 11364, 12036, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 15795, 11109, 11364, 12393, 7,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 16110, 11364, 12750, 7, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 16425, 11679, 12750, 13177, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 16866, 12036, 12750, 13653, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 17307, 12393, 12750, 14129, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 17748, 14325, 14955, 15165, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 18168, 14535, 14955, 15480, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 18588, 14745, 14955, 15795, 7,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 19008, 14955, 16110, 7, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 19428, 15165, 16110, 16425, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 20058, 15480, 16110, 16866, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 20688, 15795, 16110, 17307, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 21318, 17748, 19008, 19428, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 22018, 18168, 19008, 20058, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 22718, 18588, 19008, 20688, 7,
                                          nmax);

        simdtrf::transform_f_inner(buffer, 23418, 21318, 10, 7, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 23418, 49, nmax);

        simdtrf::transform_f_inner(buffer, 23418, 22018, 10, 7, nmax);

        simdtrf::transform_f_outer(values + 343 * nvalues + n * npairs, nvalues, buffer, 23418,
                                   49, nmax);

        simdtrf::transform_f_inner(buffer, 23418, 22718, 10, 7, nmax);

        simdtrf::transform_f_outer(values + 686 * nvalues + n * npairs, nvalues, buffer, 23418,
                                   49, nmax);
    }

    for (size_t m = 0; m < 1029; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
