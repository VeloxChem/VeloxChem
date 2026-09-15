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


#include "SimdThreeCenterElectronRepulsionGeom010RecGFD.hpp"

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
#include "SimdGeometryK1.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDF.hpp"
#include "SimdTransferDG.hpp"
#include "SimdTransferFF.hpp"
#include "SimdTransferGeom010XDF.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XDH.hpp"
#include "SimdTransferGeom010XFF.hpp"
#include "SimdTransferGeom010XFG.hpp"
#include "SimdTransferGeom010XGF.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010XPI.hpp"
#include "SimdTransferGeom010YDF.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YDH.hpp"
#include "SimdTransferGeom010YFF.hpp"
#include "SimdTransferGeom010YFG.hpp"
#include "SimdTransferGeom010YGF.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010YPI.hpp"
#include "SimdTransferGeom010ZDF.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZDH.hpp"
#include "SimdTransferGeom010ZFF.hpp"
#include "SimdTransferGeom010ZFG.hpp"
#include "SimdTransferGeom010ZGF.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferGeom010ZPI.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_gfd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_gfd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 27646, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 945 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 27646, 7267, 4264, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 10,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 18, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

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

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 48, 0, 3, 7, 8,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 54, 0, 3, 8, 9,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 60, 0, 3, 9, 10,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 66, 0, 3, 10, 11,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 72, 0, 3, 11, 12,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 78, 0, 3, 12, 13,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 84, 0, 3, 13, 14,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 90, 0, 3, 14, 15,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 96, 0, 3, 15, 16,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 102, 0, 3, 18, 21,
                                                                       48, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 112, 0, 3, 21, 24,
                                                                       54, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 122, 0, 3, 24, 27,
                                                                       60, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 132, 0, 3, 27, 30,
                                                                       66, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 142, 0, 3, 30, 33,
                                                                       72, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 152, 0, 3, 33, 36,
                                                                       78, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 162, 0, 3, 36, 39,
                                                                       84, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 172, 0, 3, 39, 42,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 182, 0, 3, 48, 54,
                                                                       102, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 197, 0, 3, 54, 60,
                                                                       112, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 212, 0, 3, 60, 66,
                                                                       122, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 227, 0, 3, 66, 72,
                                                                       132, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 242, 0, 3, 72, 78,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 257, 0, 3, 78, 84,
                                                                       152, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 272, 0, 3, 84, 90,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 287, 0, 3, 102,
                                                                       112, 182, 197, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 112,
                                                                       122, 197, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 329, 0, 3, 122,
                                                                       132, 212, 227, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 350, 0, 3, 132,
                                                                       142, 227, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 371, 0, 3, 142,
                                                                       152, 242, 257, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 392, 0, 3, 152,
                                                                       162, 257, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 413, 0, 3, 182,
                                                                       197, 287, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 441, 0, 3, 197,
                                                                       212, 308, 329, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 469, 0, 3, 212,
                                                                       227, 329, 350, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 497, 0, 3, 227,
                                                                       242, 350, 371, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 525, 0, 3, 242,
                                                                       257, 371, 392, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 553, 0, 3, 287,
                                                                       308, 413, 441, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 589, 0, 3, 308,
                                                                       329, 441, 469, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 625, 0, 3, 329,
                                                                       350, 469, 497, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 661, 0, 3, 350,
                                                                       371, 497, 525, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 697, 0, 3, 413,
                                                                       441, 553, 589, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 742, 0, 3, 441,
                                                                       469, 589, 625, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 787, 0, 3, 469,
                                                                       497, 625, 661, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 832, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 835, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 838, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 841, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 844, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 847, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 850, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 853, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 856, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 859, 3, 9, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 868, 3, 10, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 877, 3, 11, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 886, 3, 12, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 895, 3, 13, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 904, 3, 14, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 913, 3, 15, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 922, 3, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 931, 3, 24, 60,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 949, 3, 27, 66,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 967, 3, 30, 72,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 985, 3, 33, 78,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1003, 3, 36, 84,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1021, 3, 39, 90,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1039, 3, 42, 96,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1057, 3, 60, 122,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1087, 3, 66, 132,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1117, 3, 72, 142,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1147, 3, 78, 152,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1177, 3, 84, 162,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1207, 3, 90, 172,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1237, 3, 122, 212,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1282, 3, 132, 227,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1327, 3, 142, 242,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1372, 3, 152, 257,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1417, 3, 162, 272,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1462, 3, 212, 329,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1525, 3, 227, 350,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1588, 3, 242, 371,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1651, 3, 257, 392,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1714, 3, 329, 469,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1798, 3, 350, 497,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1882, 3, 371, 525,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 1966, 3, 469, 625,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 2074, 3, 497, 661,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 2182, 3, 625, 787,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2317, 3, 7, 8,
                                                                       832, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2323, 3, 8, 9,
                                                                       835, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2329, 3, 9, 10,
                                                                       838, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2335, 3, 10, 11,
                                                                       841, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2341, 3, 11, 12,
                                                                       844, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2347, 3, 12, 13,
                                                                       847, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2353, 3, 13, 14,
                                                                       850, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2359, 3, 14, 15,
                                                                       853, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2365, 3, 15, 16,
                                                                       856, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2371, 0, 3, 2317,
                                                                       832, 2323, 859, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2389, 0, 3, 2323,
                                                                       835, 2329, 868, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2407, 0, 3, 2329,
                                                                       838, 2335, 877, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2425, 0, 3, 2335,
                                                                       841, 2341, 886, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2443, 0, 3, 2341,
                                                                       844, 2347, 895, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2461, 0, 3, 2347,
                                                                       847, 2353, 904, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2479, 0, 3, 2353,
                                                                       850, 2359, 913, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2497, 0, 3, 2359,
                                                                       853, 2365, 922, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2515, 0, 3, 2371,
                                                                       859, 2389, 48, 54, 931,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2551, 0, 3, 2389,
                                                                       868, 2407, 54, 60, 949,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2587, 0, 3, 2407,
                                                                       877, 2425, 60, 66, 967,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2623, 0, 3, 2425,
                                                                       886, 2443, 66, 72, 985,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2659, 0, 3, 2443,
                                                                       895, 2461, 72, 78, 1003,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2695, 0, 3, 2461,
                                                                       904, 2479, 78, 84, 1021,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2731, 0, 3, 2479,
                                                                       913, 2497, 84, 90, 1039,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2767, 0, 3, 2515,
                                                                       931, 2551, 102, 112, 1057,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2827, 0, 3, 2551,
                                                                       949, 2587, 112, 122, 1087,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2887, 0, 3, 2587,
                                                                       967, 2623, 122, 132, 1117,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2947, 0, 3, 2623,
                                                                       985, 2659, 132, 142, 1147,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3007, 0, 3, 2659,
                                                                       1003, 2695, 142, 152,
                                                                       1177, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3067, 0, 3, 2695,
                                                                       1021, 2731, 152, 162,
                                                                       1207, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3127, 0, 3, 2767,
                                                                       1057, 2827, 182, 197,
                                                                       1237, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3217, 0, 3, 2827,
                                                                       1087, 2887, 197, 212,
                                                                       1282, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3307, 0, 3, 2887,
                                                                       1117, 2947, 212, 227,
                                                                       1327, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3397, 0, 3, 2947,
                                                                       1147, 3007, 227, 242,
                                                                       1372, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 3487, 0, 3, 3007,
                                                                       1177, 3067, 242, 257,
                                                                       1417, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3577, 0, 3, 3127,
                                                                       1237, 3217, 287, 308,
                                                                       1462, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3703, 0, 3, 3217,
                                                                       1282, 3307, 308, 329,
                                                                       1525, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3829, 0, 3, 3307,
                                                                       1327, 3397, 329, 350,
                                                                       1588, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 3955, 0, 3, 3397,
                                                                       1372, 3487, 350, 371,
                                                                       1651, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 4081, 0, 3, 3577,
                                                                       1462, 3703, 413, 441,
                                                                       1714, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 4249, 0, 3, 3703,
                                                                       1525, 3829, 441, 469,
                                                                       1798, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 4417, 0, 3, 3829,
                                                                       1588, 3955, 469, 497,
                                                                       1882, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 4585, 0, 3, 4081,
                                                                       1714, 4249, 553, 589,
                                                                       1966, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 4801, 0, 3, 4249,
                                                                       1798, 4417, 589, 625,
                                                                       2074, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 5017, 0, 3, 4585,
                                                                       1966, 4801, 697, 742,
                                                                       2182, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_f_x(buffer, 5287, 2515, 3127, 1, 6, ncols, beta);

                    simdgeo::geom_f_y(buffer, 5347, 2515, 3127, 1, 6, ncols, beta);

                    simdgeo::geom_f_z(buffer, 5407, 2515, 3127, 1, 6, ncols, beta);

                    simdgeo::geom_g_x(buffer, 5467, 2767, 3577, 1, 6, ncols, beta);

                    simdgeo::geom_g_y(buffer, 5557, 2767, 3577, 1, 6, ncols, beta);

                    simdgeo::geom_g_z(buffer, 5647, 2767, 3577, 1, 6, ncols, beta);

                    simdgeo::geom_h_x(buffer, 5737, 3127, 4081, 1, 6, ncols, beta);

                    simdgeo::geom_h_y(buffer, 5863, 3127, 4081, 1, 6, ncols, beta);

                    simdgeo::geom_h_z(buffer, 5989, 3127, 4081, 1, 6, ncols, beta);

                    simdgeo::geom_i_x(buffer, 6115, 3577, 4585, 1, 6, ncols, beta);

                    simdgeo::geom_i_y(buffer, 6283, 3577, 4585, 1, 6, ncols, beta);

                    simdgeo::geom_i_z(buffer, 6451, 3577, 4585, 1, 6, ncols, beta);

                    simdgeo::geom_k_x(buffer, 6619, 4081, 5017, 1, 6, ncols, beta);

                    simdgeo::geom_k_y(buffer, 6835, 4081, 5017, 1, 6, ncols, beta);

                    simdgeo::geom_k_z(buffer, 7051, 4081, 5017, 1, 6, ncols, beta);

                    simdfunc::contract_primitives(buffer, 7267, 5287, 60, ncols);

                    simdfunc::contract_primitives(buffer, 7377, 5347, 60, ncols);

                    simdfunc::contract_primitives(buffer, 7487, 5407, 60, ncols);

                    simdfunc::contract_primitives(buffer, 7597, 2767, 60, ncols);

                    simdfunc::contract_primitives(buffer, 7707, 5467, 90, ncols);

                    simdfunc::contract_primitives(buffer, 7872, 5557, 90, ncols);

                    simdfunc::contract_primitives(buffer, 8037, 5647, 90, ncols);

                    simdfunc::contract_primitives(buffer, 8202, 3127, 90, ncols);

                    simdfunc::contract_primitives(buffer, 8367, 5737, 126, ncols);

                    simdfunc::contract_primitives(buffer, 8598, 5863, 126, ncols);

                    simdfunc::contract_primitives(buffer, 8829, 5989, 126, ncols);

                    simdfunc::contract_primitives(buffer, 9060, 3577, 126, ncols);

                    simdfunc::contract_primitives(buffer, 9291, 6115, 168, ncols);

                    simdfunc::contract_primitives(buffer, 9599, 6283, 168, ncols);

                    simdfunc::contract_primitives(buffer, 9907, 6451, 168, ncols);

                    simdfunc::contract_primitives(buffer, 10215, 4081, 168, ncols);

                    simdfunc::contract_primitives(buffer, 10523, 6619, 216, ncols);

                    simdfunc::contract_primitives(buffer, 10919, 6835, 216, ncols);

                    simdfunc::contract_primitives(buffer, 11315, 7051, 216, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 7327, 7267, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 7437, 7377, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 7547, 7487, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 7657, 7597, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 7797, 7707, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 7962, 7872, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 8127, 8037, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 8292, 8202, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 8493, 8367, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 8724, 8598, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 8955, 8829, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 9186, 9060, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 9459, 9291, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 9767, 9599, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10075, 9907, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10383, 10215, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 10739, 10523, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 11135, 10919, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 11531, 11315, 36, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 11711, 7327, 7657, 7797, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 11861, 7437, 7657, 7962, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 12011, 7547, 7657, 8127, 5,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 12161, 7657, 8292, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 12311, 7797, 8292, 8493, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 12536, 7962, 8292, 8724, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 12761, 8127, 8292, 8955, 5,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 12986, 8292, 9186, 5, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 13211, 8493, 9186, 9459, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 13526, 8724, 9186, 9767, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 13841, 8955, 9186, 10075, 5,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 14156, 9186, 10383, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 14471, 9459, 10383, 10739, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 14891, 9767, 10383, 11135, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 15311, 10075, 10383, 11531, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 15731, 11711, 12161, 12311, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 16031, 11861, 12161, 12536, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 16331, 12011, 12161, 12761, 5,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 16631, 12161, 12986, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 16931, 12311, 12986, 13211, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 17381, 12536, 12986, 13526, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 17831, 12761, 12986, 13841, 5,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 18281, 12986, 14156, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 18731, 13211, 14156, 14471, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 19361, 13526, 14156, 14891, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 19991, 13841, 14156, 15311, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 20621, 15731, 16631, 16931, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 21121, 16031, 16631, 17381, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 21621, 16331, 16631, 17831, 5,
                                          nmax);

        simdtrf::compute_hrr_ff(buffer, coordinates, 22121, 16631, 18281, 5, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 22621, 16931, 18281, 18731, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 23371, 17381, 18281, 19361, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 24121, 17831, 18281, 19991, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gf_out_of_second(buffer, coordinates, 24871, 20621, 22121,
                                                        22621, 5, nmax);

        simdtrf::compute_hrr_geom_010y_gf_out_of_second(buffer, coordinates, 25621, 21121, 22121,
                                                        23371, 5, nmax);

        simdtrf::compute_hrr_geom_010z_gf_out_of_second(buffer, coordinates, 26371, 21621, 22121,
                                                        24121, 5, nmax);

        simdtrf::transform_f_inner(buffer, 27121, 24871, 15, 5, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 27121, 35, nmax);

        simdtrf::transform_f_inner(buffer, 27121, 25621, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 27121,
                                   35, nmax);

        simdtrf::transform_f_inner(buffer, 27121, 26371, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 630 * nvalues + n * npairs, nvalues, buffer, 27121,
                                   35, nmax);
    }

    for (size_t m = 0; m < 945; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
