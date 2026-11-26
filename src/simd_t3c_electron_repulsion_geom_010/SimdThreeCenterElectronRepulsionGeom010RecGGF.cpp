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


#include "SimdThreeCenterElectronRepulsionGeom010RecGGF.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdGeometryG1.hpp"
#include "SimdGeometryH1.hpp"
#include "SimdGeometryI1.hpp"
#include "SimdGeometryK1.hpp"
#include "SimdGeometryL1.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDG.hpp"
#include "SimdTransferDH.hpp"
#include "SimdTransferFG.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XDH.hpp"
#include "SimdTransferGeom010XDI.hpp"
#include "SimdTransferGeom010XFG.hpp"
#include "SimdTransferGeom010XFH.hpp"
#include "SimdTransferGeom010XGG.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010XPI.hpp"
#include "SimdTransferGeom010XPK.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YDH.hpp"
#include "SimdTransferGeom010YDI.hpp"
#include "SimdTransferGeom010YFG.hpp"
#include "SimdTransferGeom010YFH.hpp"
#include "SimdTransferGeom010YGG.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010YPI.hpp"
#include "SimdTransferGeom010YPK.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZDH.hpp"
#include "SimdTransferGeom010ZDI.hpp"
#include "SimdTransferGeom010ZFG.hpp"
#include "SimdTransferGeom010ZFH.hpp"
#include "SimdTransferGeom010ZGG.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferGeom010ZPI.hpp"
#include "SimdTransferGeom010ZPK.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_ggf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_ggf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 60994, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1701 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 60994, 20399, 8780, dimensions);

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
                                                        5, 6, 7, 8, 9, 10, 11, 12}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 19, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 52, 0, 3, 7, 8,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 58, 0, 3, 8, 9,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 64, 0, 3, 9, 10,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 70, 0, 3, 10, 11,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 76, 0, 3, 11, 12,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 82, 0, 3, 12, 13,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 88, 0, 3, 13, 14,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 94, 0, 3, 14, 15,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 100, 0, 3, 15, 16,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 106, 0, 3, 16, 17,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 112, 0, 3, 19, 22,
                                                                       52, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 122, 0, 3, 22, 25,
                                                                       58, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 132, 0, 3, 25, 28,
                                                                       64, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 142, 0, 3, 28, 31,
                                                                       70, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 152, 0, 3, 31, 34,
                                                                       76, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 162, 0, 3, 34, 37,
                                                                       82, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 172, 0, 3, 37, 40,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 182, 0, 3, 40, 43,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 192, 0, 3, 43, 46,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 202, 0, 3, 52, 58,
                                                                       112, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 217, 0, 3, 58, 64,
                                                                       122, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 232, 0, 3, 64, 70,
                                                                       132, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 247, 0, 3, 70, 76,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 262, 0, 3, 76, 82,
                                                                       152, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 277, 0, 3, 82, 88,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 292, 0, 3, 88, 94,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 307, 0, 3, 94,
                                                                       100, 182, 192, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 322, 0, 3, 112,
                                                                       122, 202, 217, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 343, 0, 3, 122,
                                                                       132, 217, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 364, 0, 3, 132,
                                                                       142, 232, 247, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 385, 0, 3, 142,
                                                                       152, 247, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 406, 0, 3, 152,
                                                                       162, 262, 277, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 427, 0, 3, 162,
                                                                       172, 277, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 172,
                                                                       182, 292, 307, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 469, 0, 3, 202,
                                                                       217, 322, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 497, 0, 3, 217,
                                                                       232, 343, 364, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 525, 0, 3, 232,
                                                                       247, 364, 385, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 553, 0, 3, 247,
                                                                       262, 385, 406, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 581, 0, 3, 262,
                                                                       277, 406, 427, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 609, 0, 3, 277,
                                                                       292, 427, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 637, 0, 3, 322,
                                                                       343, 469, 497, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 673, 0, 3, 343,
                                                                       364, 497, 525, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 709, 0, 3, 364,
                                                                       385, 525, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 745, 0, 3, 385,
                                                                       406, 553, 581, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 781, 0, 3, 406,
                                                                       427, 581, 609, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 817, 0, 3, 469,
                                                                       497, 637, 673, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 862, 0, 3, 497,
                                                                       525, 673, 709, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 907, 0, 3, 525,
                                                                       553, 709, 745, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 952, 0, 3, 553,
                                                                       581, 745, 781, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 997, 0, 3, 637,
                                                                       673, 817, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1052, 0, 3, 673,
                                                                       709, 862, 907, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1107, 0, 3, 709,
                                                                       745, 907, 952, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1162, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1165, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1168, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1171, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1174, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1177, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1180, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1183, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1186, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1189, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1192, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1195, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1198, 3, 9, 25,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1207, 3, 10, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1216, 3, 11, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1225, 3, 12, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1234, 3, 13, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1243, 3, 14, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1252, 3, 15, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1261, 3, 16, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1270, 3, 17, 49,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1279, 3, 19, 52,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1297, 3, 22, 58,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1315, 3, 25, 64,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1333, 3, 28, 70,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1351, 3, 31, 76,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1369, 3, 34, 82,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1387, 3, 37, 88,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1405, 3, 40, 94,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1423, 3, 43, 100,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1441, 3, 46, 106,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1459, 3, 52, 112,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1489, 3, 58, 122,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1519, 3, 64, 132,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1549, 3, 70, 142,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1579, 3, 76, 152,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1609, 3, 82, 162,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1639, 3, 88, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1669, 3, 94, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1699, 3, 100, 192,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1729, 3, 112, 202,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1774, 3, 122, 217,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1819, 3, 132, 232,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1864, 3, 142, 247,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1909, 3, 152, 262,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1954, 3, 162, 277,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1999, 3, 172, 292,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2044, 3, 182, 307,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2089, 3, 202, 322,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2152, 3, 217, 343,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2215, 3, 232, 364,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2278, 3, 247, 385,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2341, 3, 262, 406,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2404, 3, 277, 427,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2467, 3, 292, 448,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2530, 3, 322, 469,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2614, 3, 343, 497,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2698, 3, 364, 525,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2782, 3, 385, 553,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2866, 3, 406, 581,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2950, 3, 427, 609,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3034, 3, 469, 637,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3142, 3, 497, 673,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3250, 3, 525, 709,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3358, 3, 553, 745,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3466, 3, 581, 781,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3574, 3, 637, 817,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3709, 3, 673, 862,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3844, 3, 709, 907,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3979, 3, 745, 952,
                                                                       ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 4114, 3, 817, 997,
                                                                       ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 4279, 3, 862,
                                                                       1052, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 4444, 3, 907,
                                                                       1107, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4609, 3, 7, 8,
                                                                       1168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4615, 3, 8, 9,
                                                                       1171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4621, 3, 9, 10,
                                                                       1174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4627, 3, 10, 11,
                                                                       1177, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4633, 3, 11, 12,
                                                                       1180, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4639, 3, 12, 13,
                                                                       1183, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4645, 3, 13, 14,
                                                                       1186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4651, 3, 14, 15,
                                                                       1189, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4657, 3, 15, 16,
                                                                       1192, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4663, 3, 16, 17,
                                                                       1195, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4669, 0, 3, 4609,
                                                                       1168, 4615, 1198, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4687, 0, 3, 4615,
                                                                       1171, 4621, 1207, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4705, 0, 3, 4621,
                                                                       1174, 4627, 1216, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4723, 0, 3, 4627,
                                                                       1177, 4633, 1225, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4741, 0, 3, 4633,
                                                                       1180, 4639, 1234, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4759, 0, 3, 4639,
                                                                       1183, 4645, 1243, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4777, 0, 3, 4645,
                                                                       1186, 4651, 1252, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4795, 0, 3, 4651,
                                                                       1189, 4657, 1261, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4813, 0, 3, 4657,
                                                                       1192, 4663, 1270, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4831, 0, 3, 4669,
                                                                       1198, 4687, 52, 58, 1315,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4867, 0, 3, 4687,
                                                                       1207, 4705, 58, 64, 1333,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4903, 0, 3, 4705,
                                                                       1216, 4723, 64, 70, 1351,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4939, 0, 3, 4723,
                                                                       1225, 4741, 70, 76, 1369,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4975, 0, 3, 4741,
                                                                       1234, 4759, 76, 82, 1387,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5011, 0, 3, 4759,
                                                                       1243, 4777, 82, 88, 1405,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5047, 0, 3, 4777,
                                                                       1252, 4795, 88, 94, 1423,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5083, 0, 3, 4795,
                                                                       1261, 4813, 94, 100, 1441,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5119, 0, 3, 4831,
                                                                       1315, 4867, 112, 122,
                                                                       1519, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5179, 0, 3, 4867,
                                                                       1333, 4903, 122, 132,
                                                                       1549, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5239, 0, 3, 4903,
                                                                       1351, 4939, 132, 142,
                                                                       1579, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5299, 0, 3, 4939,
                                                                       1369, 4975, 142, 152,
                                                                       1609, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5359, 0, 3, 4975,
                                                                       1387, 5011, 152, 162,
                                                                       1639, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5419, 0, 3, 5011,
                                                                       1405, 5047, 162, 172,
                                                                       1669, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5479, 0, 3, 5047,
                                                                       1423, 5083, 172, 182,
                                                                       1699, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5539, 0, 3, 5119,
                                                                       1519, 5179, 202, 217,
                                                                       1819, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5629, 0, 3, 5179,
                                                                       1549, 5239, 217, 232,
                                                                       1864, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5719, 0, 3, 5239,
                                                                       1579, 5299, 232, 247,
                                                                       1909, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5809, 0, 3, 5299,
                                                                       1609, 5359, 247, 262,
                                                                       1954, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5899, 0, 3, 5359,
                                                                       1639, 5419, 262, 277,
                                                                       1999, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5989, 0, 3, 5419,
                                                                       1669, 5479, 277, 292,
                                                                       2044, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6079, 0, 3, 5539,
                                                                       1819, 5629, 322, 343,
                                                                       2215, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6205, 0, 3, 5629,
                                                                       1864, 5719, 343, 364,
                                                                       2278, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6331, 0, 3, 5719,
                                                                       1909, 5809, 364, 385,
                                                                       2341, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6457, 0, 3, 5809,
                                                                       1954, 5899, 385, 406,
                                                                       2404, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6583, 0, 3, 5899,
                                                                       1999, 5989, 406, 427,
                                                                       2467, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6709, 0, 3, 6079,
                                                                       2215, 6205, 469, 497,
                                                                       2698, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6877, 0, 3, 6205,
                                                                       2278, 6331, 497, 525,
                                                                       2782, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7045, 0, 3, 6331,
                                                                       2341, 6457, 525, 553,
                                                                       2866, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7213, 0, 3, 6457,
                                                                       2404, 6583, 553, 581,
                                                                       2950, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 7381, 0, 3, 6709,
                                                                       2698, 6877, 637, 673,
                                                                       3250, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 7597, 0, 3, 6877,
                                                                       2782, 7045, 673, 709,
                                                                       3358, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 7813, 0, 3, 7045,
                                                                       2866, 7213, 709, 745,
                                                                       3466, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 8029, 0, 3, 7381,
                                                                       3250, 7597, 817, 862,
                                                                       3844, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 8299, 0, 3, 7597,
                                                                       3358, 7813, 862, 907,
                                                                       3979, ncols, gamma, p,
                                                                       q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 8569, 0, 3, 8029,
                                                                       3844, 8299, 997, 1052,
                                                                       4444, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8899, 3, 1162,
                                                                       1165, 4609, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8909, 3, 1165,
                                                                       1168, 4615, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8919, 3, 1168,
                                                                       1171, 4621, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8929, 3, 1171,
                                                                       1174, 4627, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8939, 3, 1174,
                                                                       1177, 4633, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8949, 3, 1177,
                                                                       1180, 4639, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8959, 3, 1180,
                                                                       1183, 4645, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8969, 3, 1183,
                                                                       1186, 4651, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8979, 3, 1186,
                                                                       1189, 4657, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8989, 3, 1189,
                                                                       1192, 4663, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8999, 0, 3, 8899,
                                                                       4609, 8909, 4669, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9029, 0, 3, 8909,
                                                                       4615, 8919, 4687, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9059, 0, 3, 8919,
                                                                       4621, 8929, 4705, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9089, 0, 3, 8929,
                                                                       4627, 8939, 4723, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9119, 0, 3, 8939,
                                                                       4633, 8949, 4741, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9149, 0, 3, 8949,
                                                                       4639, 8959, 4759, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9179, 0, 3, 8959,
                                                                       4645, 8969, 4777, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9209, 0, 3, 8969,
                                                                       4651, 8979, 4795, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9239, 0, 3, 8979,
                                                                       4657, 8989, 4813, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9269, 0, 3, 8999,
                                                                       4669, 9029, 1279, 1297,
                                                                       4831, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9329, 0, 3, 9029,
                                                                       4687, 9059, 1297, 1315,
                                                                       4867, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9389, 0, 3, 9059,
                                                                       4705, 9089, 1315, 1333,
                                                                       4903, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9449, 0, 3, 9089,
                                                                       4723, 9119, 1333, 1351,
                                                                       4939, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9509, 0, 3, 9119,
                                                                       4741, 9149, 1351, 1369,
                                                                       4975, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9569, 0, 3, 9149,
                                                                       4759, 9179, 1369, 1387,
                                                                       5011, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9629, 0, 3, 9179,
                                                                       4777, 9209, 1387, 1405,
                                                                       5047, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9689, 0, 3, 9209,
                                                                       4795, 9239, 1405, 1423,
                                                                       5083, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9749, 0, 3, 9269,
                                                                       4831, 9329, 1459, 1489,
                                                                       5119, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9849, 0, 3, 9329,
                                                                       4867, 9389, 1489, 1519,
                                                                       5179, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9949, 0, 3, 9389,
                                                                       4903, 9449, 1519, 1549,
                                                                       5239, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10049, 0, 3, 9449,
                                                                       4939, 9509, 1549, 1579,
                                                                       5299, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10149, 0, 3, 9509,
                                                                       4975, 9569, 1579, 1609,
                                                                       5359, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10249, 0, 3, 9569,
                                                                       5011, 9629, 1609, 1639,
                                                                       5419, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10349, 0, 3, 9629,
                                                                       5047, 9689, 1639, 1669,
                                                                       5479, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10449, 0, 3, 9749,
                                                                       5119, 9849, 1729, 1774,
                                                                       5539, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10599, 0, 3, 9849,
                                                                       5179, 9949, 1774, 1819,
                                                                       5629, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10749, 0, 3, 9949,
                                                                       5239, 10049, 1819, 1864,
                                                                       5719, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10899, 0, 3,
                                                                       10049, 5299, 10149, 1864,
                                                                       1909, 5809, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11049, 0, 3,
                                                                       10149, 5359, 10249, 1909,
                                                                       1954, 5899, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11199, 0, 3,
                                                                       10249, 5419, 10349, 1954,
                                                                       1999, 5989, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 11349, 0, 3,
                                                                       10449, 5539, 10599, 2089,
                                                                       2152, 6079, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 11559, 0, 3,
                                                                       10599, 5629, 10749, 2152,
                                                                       2215, 6205, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 11769, 0, 3,
                                                                       10749, 5719, 10899, 2215,
                                                                       2278, 6331, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 11979, 0, 3,
                                                                       10899, 5809, 11049, 2278,
                                                                       2341, 6457, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12189, 0, 3,
                                                                       11049, 5899, 11199, 2341,
                                                                       2404, 6583, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 12399, 0, 3,
                                                                       11349, 6079, 11559, 2530,
                                                                       2614, 6709, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 12679, 0, 3,
                                                                       11559, 6205, 11769, 2614,
                                                                       2698, 6877, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 12959, 0, 3,
                                                                       11769, 6331, 11979, 2698,
                                                                       2782, 7045, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 13239, 0, 3,
                                                                       11979, 6457, 12189, 2782,
                                                                       2866, 7213, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 13519, 0, 3,
                                                                       12399, 6709, 12679, 3034,
                                                                       3142, 7381, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 13879, 0, 3,
                                                                       12679, 6877, 12959, 3142,
                                                                       3250, 7597, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 14239, 0, 3,
                                                                       12959, 7045, 13239, 3250,
                                                                       3358, 7813, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 14599, 0, 3,
                                                                       13519, 7381, 13879, 3574,
                                                                       3709, 8029, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 15049, 0, 3,
                                                                       13879, 7597, 14239, 3709,
                                                                       3844, 8299, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 15499, 0, 3,
                                                                       14599, 8029, 15049, 4114,
                                                                       4279, 8569, ncols, gamma,
                                                                       p, q);

                    simdgeo::geom_g_x(buffer, 16049, 9749, 11349, 1, 10, ncols, beta);

                    simdgeo::geom_g_y(buffer, 16199, 9749, 11349, 1, 10, ncols, beta);

                    simdgeo::geom_g_z(buffer, 16349, 9749, 11349, 1, 10, ncols, beta);

                    simdgeo::geom_h_x(buffer, 16499, 10449, 12399, 1, 10, ncols, beta);

                    simdgeo::geom_h_y(buffer, 16709, 10449, 12399, 1, 10, ncols, beta);

                    simdgeo::geom_h_z(buffer, 16919, 10449, 12399, 1, 10, ncols, beta);

                    simdgeo::geom_i_x(buffer, 17129, 11349, 13519, 1, 10, ncols, beta);

                    simdgeo::geom_i_y(buffer, 17409, 11349, 13519, 1, 10, ncols, beta);

                    simdgeo::geom_i_z(buffer, 17689, 11349, 13519, 1, 10, ncols, beta);

                    simdgeo::geom_k_x(buffer, 17969, 12399, 14599, 1, 10, ncols, beta);

                    simdgeo::geom_k_y(buffer, 18329, 12399, 14599, 1, 10, ncols, beta);

                    simdgeo::geom_k_z(buffer, 18689, 12399, 14599, 1, 10, ncols, beta);

                    simdgeo::geom_l_x(buffer, 19049, 13519, 15499, 1, 10, ncols, beta);

                    simdgeo::geom_l_y(buffer, 19499, 13519, 15499, 1, 10, ncols, beta);

                    simdgeo::geom_l_z(buffer, 19949, 13519, 15499, 1, 10, ncols, beta);

                    simdfunc::contract_primitives(buffer, 20399, 16049, 150, ncols);

                    simdfunc::contract_primitives(buffer, 20654, 16199, 150, ncols);

                    simdfunc::contract_primitives(buffer, 20909, 16349, 150, ncols);

                    simdfunc::contract_primitives(buffer, 21164, 10449, 150, ncols);

                    simdfunc::contract_primitives(buffer, 21419, 16499, 210, ncols);

                    simdfunc::contract_primitives(buffer, 21776, 16709, 210, ncols);

                    simdfunc::contract_primitives(buffer, 22133, 16919, 210, ncols);

                    simdfunc::contract_primitives(buffer, 22490, 11349, 210, ncols);

                    simdfunc::contract_primitives(buffer, 22847, 17129, 280, ncols);

                    simdfunc::contract_primitives(buffer, 23323, 17409, 280, ncols);

                    simdfunc::contract_primitives(buffer, 23799, 17689, 280, ncols);

                    simdfunc::contract_primitives(buffer, 24275, 12399, 280, ncols);

                    simdfunc::contract_primitives(buffer, 24751, 17969, 360, ncols);

                    simdfunc::contract_primitives(buffer, 25363, 18329, 360, ncols);

                    simdfunc::contract_primitives(buffer, 25975, 18689, 360, ncols);

                    simdfunc::contract_primitives(buffer, 26587, 13519, 360, ncols);

                    simdfunc::contract_primitives(buffer, 27199, 19049, 450, ncols);

                    simdfunc::contract_primitives(buffer, 27964, 19499, 450, ncols);

                    simdfunc::contract_primitives(buffer, 28729, 19949, 450, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 20549, 20399, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 20804, 20654, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 21059, 20909, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 21314, 21164, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 21629, 21419, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 21986, 21776, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 22343, 22133, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 22700, 22490, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 23127, 22847, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 23603, 23323, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 24079, 23799, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 24555, 24275, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 25111, 24751, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 25723, 25363, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 26335, 25975, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 26947, 26587, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 27649, 27199, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 28414, 27964, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 29179, 28729, 45, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 29494, 20549, 21314, 21629, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 29809, 20804, 21314, 21986, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 30124, 21059, 21314, 22343, 7,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 30439, 21314, 22700, 7, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 30754, 21629, 22700, 23127, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 31195, 21986, 22700, 23603, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 31636, 22343, 22700, 24079, 7,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 32077, 22700, 24555, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 32518, 23127, 24555, 25111, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 33106, 23603, 24555, 25723, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 33694, 24079, 24555, 26335, 7,
                                          nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 34282, 24555, 26947, 7, nmax);

        simdtrf::compute_hrr_geom_010x_pk(buffer, coordinates, 34870, 25111, 26947, 27649, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pk(buffer, coordinates, 35626, 25723, 26947, 28414, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pk(buffer, coordinates, 36382, 26335, 26947, 29179, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 37138, 29494, 30439, 30754, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 37768, 29809, 30439, 31195, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 38398, 30124, 30439, 31636, 7,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 39028, 30439, 32077, 7, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 39658, 30754, 32077, 32518, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 40540, 31195, 32077, 33106, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 41422, 31636, 32077, 33694, 7,
                                          nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 42304, 32077, 34282, 7, nmax);

        simdtrf::compute_hrr_geom_010x_di(buffer, coordinates, 43186, 32518, 34282, 34870, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_di(buffer, coordinates, 44362, 33106, 34282, 35626, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_di(buffer, coordinates, 45538, 33694, 34282, 36382, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 46714, 37138, 39028, 39658, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 47764, 37768, 39028, 40540, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 48814, 38398, 39028, 41422, 7,
                                          nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 49864, 39028, 42304, 7, nmax);

        simdtrf::compute_hrr_geom_010x_fh(buffer, coordinates, 50914, 39658, 42304, 43186, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fh(buffer, coordinates, 52384, 40540, 42304, 44362, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fh(buffer, coordinates, 53854, 41422, 42304, 45538, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gg(buffer, coordinates, 55324, 46714, 49864, 50914, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_gg(buffer, coordinates, 56899, 47764, 49864, 52384, 7,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_gg(buffer, coordinates, 58474, 48814, 49864, 53854, 7,
                                          nmax);

        simdtrf::transform_g_inner(buffer, 60049, 55324, 15, 7, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 60049, 63, nmax);

        simdtrf::transform_g_inner(buffer, 60049, 56899, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 567 * nvalues + n * npairs, nvalues, buffer, 60049,
                                   63, nmax);

        simdtrf::transform_g_inner(buffer, 60049, 58474, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 1134 * nvalues + n * npairs, nvalues, buffer, 60049,
                                   63, nmax);
    }

    for (size_t m = 0; m < 1701; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
