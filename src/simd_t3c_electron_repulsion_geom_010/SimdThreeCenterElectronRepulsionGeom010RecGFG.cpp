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


#include "SimdThreeCenterElectronRepulsionGeom010RecGFG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
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
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_gfg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_gfg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 64291, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 64291, 25912, 9372, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 12,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 20, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 56, 0, 3, 7, 8,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 62, 0, 3, 8, 9,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 68, 0, 3, 9, 10,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 74, 0, 3, 10, 11,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 80, 0, 3, 11, 12,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 86, 0, 3, 12, 13,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 92, 0, 3, 13, 14,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 98, 0, 3, 14, 15,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 104, 0, 3, 15, 16,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 110, 0, 3, 16, 17,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 116, 0, 3, 17, 18,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 122, 0, 3, 20, 23,
                                                                       56, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 132, 0, 3, 23, 26,
                                                                       62, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 142, 0, 3, 26, 29,
                                                                       68, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 152, 0, 3, 29, 32,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 162, 0, 3, 32, 35,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 172, 0, 3, 35, 38,
                                                                       86, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 182, 0, 3, 38, 41,
                                                                       92, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 192, 0, 3, 41, 44,
                                                                       98, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 202, 0, 3, 44, 47,
                                                                       104, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 212, 0, 3, 47, 50,
                                                                       110, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 222, 0, 3, 56, 62,
                                                                       122, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 237, 0, 3, 62, 68,
                                                                       132, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 252, 0, 3, 68, 74,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 267, 0, 3, 74, 80,
                                                                       152, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 282, 0, 3, 80, 86,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 297, 0, 3, 86, 92,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 312, 0, 3, 92, 98,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 327, 0, 3, 98,
                                                                       104, 192, 202, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 342, 0, 3, 104,
                                                                       110, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 357, 0, 3, 122,
                                                                       132, 222, 237, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 132,
                                                                       142, 237, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 399, 0, 3, 142,
                                                                       152, 252, 267, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 420, 0, 3, 152,
                                                                       162, 267, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 441, 0, 3, 162,
                                                                       172, 282, 297, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 462, 0, 3, 172,
                                                                       182, 297, 312, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 483, 0, 3, 182,
                                                                       192, 312, 327, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 504, 0, 3, 192,
                                                                       202, 327, 342, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 525, 0, 3, 222,
                                                                       237, 357, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 553, 0, 3, 237,
                                                                       252, 378, 399, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 581, 0, 3, 252,
                                                                       267, 399, 420, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 609, 0, 3, 267,
                                                                       282, 420, 441, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 637, 0, 3, 282,
                                                                       297, 441, 462, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 665, 0, 3, 297,
                                                                       312, 462, 483, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 693, 0, 3, 312,
                                                                       327, 483, 504, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 721, 0, 3, 357,
                                                                       378, 525, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 757, 0, 3, 378,
                                                                       399, 553, 581, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 793, 0, 3, 399,
                                                                       420, 581, 609, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 829, 0, 3, 420,
                                                                       441, 609, 637, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 865, 0, 3, 441,
                                                                       462, 637, 665, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 901, 0, 3, 462,
                                                                       483, 665, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 937, 0, 3, 525,
                                                                       553, 721, 757, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 982, 0, 3, 553,
                                                                       581, 757, 793, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1027, 0, 3, 581,
                                                                       609, 793, 829, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 609,
                                                                       637, 829, 865, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1117, 0, 3, 637,
                                                                       665, 865, 901, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1162, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1165, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1168, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1171, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1174, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1177, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1180, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1183, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1186, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1189, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1192, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1195, 3, 9, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1204, 3, 10, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1213, 3, 11, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1222, 3, 12, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1231, 3, 13, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1240, 3, 14, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1249, 3, 15, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1258, 3, 16, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1267, 3, 17, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1276, 3, 18, 53,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1285, 3, 26, 68,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1303, 3, 29, 74,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1321, 3, 32, 80,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1339, 3, 35, 86,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1357, 3, 38, 92,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1375, 3, 41, 98,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1393, 3, 44, 104,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1411, 3, 47, 110,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1429, 3, 50, 116,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1447, 3, 68, 142,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1477, 3, 74, 152,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1507, 3, 80, 162,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1537, 3, 86, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1567, 3, 92, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1597, 3, 98, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1627, 3, 104, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1657, 3, 110, 212,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1687, 3, 142, 252,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1732, 3, 152, 267,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1777, 3, 162, 282,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1822, 3, 172, 297,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1867, 3, 182, 312,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1912, 3, 192, 327,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1957, 3, 202, 342,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2002, 3, 252, 399,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2065, 3, 267, 420,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2128, 3, 282, 441,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2191, 3, 297, 462,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2254, 3, 312, 483,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2317, 3, 327, 504,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2380, 3, 399, 581,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2464, 3, 420, 609,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2548, 3, 441, 637,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2632, 3, 462, 665,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2716, 3, 483, 693,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 2800, 3, 581, 793,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 2908, 3, 609, 829,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3016, 3, 637, 865,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3124, 3, 665, 901,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3232, 3, 793,
                                                                       1027, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3367, 3, 829,
                                                                       1072, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 3502, 3, 865,
                                                                       1117, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3637, 3, 7, 8,
                                                                       1162, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3643, 3, 8, 9,
                                                                       1165, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3649, 3, 9, 10,
                                                                       1168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3655, 3, 10, 11,
                                                                       1171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3661, 3, 11, 12,
                                                                       1174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3667, 3, 12, 13,
                                                                       1177, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3673, 3, 13, 14,
                                                                       1180, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3679, 3, 14, 15,
                                                                       1183, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3685, 3, 15, 16,
                                                                       1186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3691, 3, 16, 17,
                                                                       1189, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3697, 3, 17, 18,
                                                                       1192, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3703, 0, 3, 3637,
                                                                       1162, 3643, 1195, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3721, 0, 3, 3643,
                                                                       1165, 3649, 1204, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3739, 0, 3, 3649,
                                                                       1168, 3655, 1213, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3757, 0, 3, 3655,
                                                                       1171, 3661, 1222, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3775, 0, 3, 3661,
                                                                       1174, 3667, 1231, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3793, 0, 3, 3667,
                                                                       1177, 3673, 1240, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3811, 0, 3, 3673,
                                                                       1180, 3679, 1249, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3829, 0, 3, 3679,
                                                                       1183, 3685, 1258, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3847, 0, 3, 3685,
                                                                       1186, 3691, 1267, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3865, 0, 3, 3691,
                                                                       1189, 3697, 1276, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3883, 0, 3, 3703,
                                                                       1195, 3721, 56, 62, 1285,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3919, 0, 3, 3721,
                                                                       1204, 3739, 62, 68, 1303,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3955, 0, 3, 3739,
                                                                       1213, 3757, 68, 74, 1321,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3991, 0, 3, 3757,
                                                                       1222, 3775, 74, 80, 1339,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4027, 0, 3, 3775,
                                                                       1231, 3793, 80, 86, 1357,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4063, 0, 3, 3793,
                                                                       1240, 3811, 86, 92, 1375,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4099, 0, 3, 3811,
                                                                       1249, 3829, 92, 98, 1393,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4135, 0, 3, 3829,
                                                                       1258, 3847, 98, 104, 1411,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4171, 0, 3, 3847,
                                                                       1267, 3865, 104, 110,
                                                                       1429, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4207, 0, 3, 3883,
                                                                       1285, 3919, 122, 132,
                                                                       1447, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4267, 0, 3, 3919,
                                                                       1303, 3955, 132, 142,
                                                                       1477, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4327, 0, 3, 3955,
                                                                       1321, 3991, 142, 152,
                                                                       1507, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4387, 0, 3, 3991,
                                                                       1339, 4027, 152, 162,
                                                                       1537, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4447, 0, 3, 4027,
                                                                       1357, 4063, 162, 172,
                                                                       1567, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4507, 0, 3, 4063,
                                                                       1375, 4099, 172, 182,
                                                                       1597, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4567, 0, 3, 4099,
                                                                       1393, 4135, 182, 192,
                                                                       1627, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4627, 0, 3, 4135,
                                                                       1411, 4171, 192, 202,
                                                                       1657, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4687, 0, 3, 4207,
                                                                       1447, 4267, 222, 237,
                                                                       1687, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4777, 0, 3, 4267,
                                                                       1477, 4327, 237, 252,
                                                                       1732, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4867, 0, 3, 4327,
                                                                       1507, 4387, 252, 267,
                                                                       1777, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4957, 0, 3, 4387,
                                                                       1537, 4447, 267, 282,
                                                                       1822, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5047, 0, 3, 4447,
                                                                       1567, 4507, 282, 297,
                                                                       1867, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5137, 0, 3, 4507,
                                                                       1597, 4567, 297, 312,
                                                                       1912, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5227, 0, 3, 4567,
                                                                       1627, 4627, 312, 327,
                                                                       1957, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5317, 0, 3, 4687,
                                                                       1687, 4777, 357, 378,
                                                                       2002, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5443, 0, 3, 4777,
                                                                       1732, 4867, 378, 399,
                                                                       2065, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5569, 0, 3, 4867,
                                                                       1777, 4957, 399, 420,
                                                                       2128, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5695, 0, 3, 4957,
                                                                       1822, 5047, 420, 441,
                                                                       2191, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5821, 0, 3, 5047,
                                                                       1867, 5137, 441, 462,
                                                                       2254, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 5947, 0, 3, 5137,
                                                                       1912, 5227, 462, 483,
                                                                       2317, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6073, 0, 3, 5317,
                                                                       2002, 5443, 525, 553,
                                                                       2380, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6241, 0, 3, 5443,
                                                                       2065, 5569, 553, 581,
                                                                       2464, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6409, 0, 3, 5569,
                                                                       2128, 5695, 581, 609,
                                                                       2548, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6577, 0, 3, 5695,
                                                                       2191, 5821, 609, 637,
                                                                       2632, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6745, 0, 3, 5821,
                                                                       2254, 5947, 637, 665,
                                                                       2716, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 6913, 0, 3, 6073,
                                                                       2380, 6241, 721, 757,
                                                                       2800, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 7129, 0, 3, 6241,
                                                                       2464, 6409, 757, 793,
                                                                       2908, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 7345, 0, 3, 6409,
                                                                       2548, 6577, 793, 829,
                                                                       3016, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 7561, 0, 3, 6577,
                                                                       2632, 6745, 829, 865,
                                                                       3124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 7777, 0, 3, 6913,
                                                                       2800, 7129, 937, 982,
                                                                       3232, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 8047, 0, 3, 7129,
                                                                       2908, 7345, 982, 1027,
                                                                       3367, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 8317, 0, 3, 7345,
                                                                       3016, 7561, 1027, 1072,
                                                                       3502, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8587, 3, 1162,
                                                                       1165, 3649, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8597, 3, 1165,
                                                                       1168, 3655, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8607, 3, 1168,
                                                                       1171, 3661, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8617, 3, 1171,
                                                                       1174, 3667, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8627, 3, 1174,
                                                                       1177, 3673, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8637, 3, 1177,
                                                                       1180, 3679, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8647, 3, 1180,
                                                                       1183, 3685, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8657, 3, 1183,
                                                                       1186, 3691, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8667, 3, 1186,
                                                                       1189, 3697, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8677, 0, 3, 8587,
                                                                       3649, 8597, 3739, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8707, 0, 3, 8597,
                                                                       3655, 8607, 3757, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8737, 0, 3, 8607,
                                                                       3661, 8617, 3775, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8767, 0, 3, 8617,
                                                                       3667, 8627, 3793, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8797, 0, 3, 8627,
                                                                       3673, 8637, 3811, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8827, 0, 3, 8637,
                                                                       3679, 8647, 3829, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8857, 0, 3, 8647,
                                                                       3685, 8657, 3847, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8887, 0, 3, 8657,
                                                                       3691, 8667, 3865, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8917, 0, 3, 8677,
                                                                       3739, 8707, 1285, 1303,
                                                                       3955, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8977, 0, 3, 8707,
                                                                       3757, 8737, 1303, 1321,
                                                                       3991, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9037, 0, 3, 8737,
                                                                       3775, 8767, 1321, 1339,
                                                                       4027, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9097, 0, 3, 8767,
                                                                       3793, 8797, 1339, 1357,
                                                                       4063, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9157, 0, 3, 8797,
                                                                       3811, 8827, 1357, 1375,
                                                                       4099, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9217, 0, 3, 8827,
                                                                       3829, 8857, 1375, 1393,
                                                                       4135, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9277, 0, 3, 8857,
                                                                       3847, 8887, 1393, 1411,
                                                                       4171, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9337, 0, 3, 8917,
                                                                       3955, 8977, 1447, 1477,
                                                                       4327, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9437, 0, 3, 8977,
                                                                       3991, 9037, 1477, 1507,
                                                                       4387, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9537, 0, 3, 9037,
                                                                       4027, 9097, 1507, 1537,
                                                                       4447, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9637, 0, 3, 9097,
                                                                       4063, 9157, 1537, 1567,
                                                                       4507, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9737, 0, 3, 9157,
                                                                       4099, 9217, 1567, 1597,
                                                                       4567, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9837, 0, 3, 9217,
                                                                       4135, 9277, 1597, 1627,
                                                                       4627, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 9937, 0, 3, 9337,
                                                                       4327, 9437, 1687, 1732,
                                                                       4867, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10087, 0, 3, 9437,
                                                                       4387, 9537, 1732, 1777,
                                                                       4957, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10237, 0, 3, 9537,
                                                                       4447, 9637, 1777, 1822,
                                                                       5047, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10387, 0, 3, 9637,
                                                                       4507, 9737, 1822, 1867,
                                                                       5137, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10537, 0, 3, 9737,
                                                                       4567, 9837, 1867, 1912,
                                                                       5227, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 10687, 0, 3, 9937,
                                                                       4867, 10087, 2002, 2065,
                                                                       5569, ncols, gamma, p,
                                                                       q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 10897, 0, 3,
                                                                       10087, 4957, 10237, 2065,
                                                                       2128, 5695, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 11107, 0, 3,
                                                                       10237, 5047, 10387, 2128,
                                                                       2191, 5821, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 11317, 0, 3,
                                                                       10387, 5137, 10537, 2191,
                                                                       2254, 5947, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 11527, 0, 3,
                                                                       10687, 5569, 10897, 2380,
                                                                       2464, 6409, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 11807, 0, 3,
                                                                       10897, 5695, 11107, 2464,
                                                                       2548, 6577, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 12087, 0, 3,
                                                                       11107, 5821, 11317, 2548,
                                                                       2632, 6745, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 12367, 0, 3,
                                                                       11527, 6409, 11807, 2800,
                                                                       2908, 7345, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 12727, 0, 3,
                                                                       11807, 6577, 12087, 2908,
                                                                       3016, 7561, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 13087, 0, 3,
                                                                       12367, 7345, 12727, 3232,
                                                                       3367, 8317, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13537, 3, 3637,
                                                                       3643, 8587, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13552, 3, 3643,
                                                                       3649, 8597, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13567, 3, 3649,
                                                                       3655, 8607, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13582, 3, 3655,
                                                                       3661, 8617, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13597, 3, 3661,
                                                                       3667, 8627, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13612, 3, 3667,
                                                                       3673, 8637, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13627, 3, 3673,
                                                                       3679, 8647, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13642, 3, 3679,
                                                                       3685, 8657, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13657, 3, 3685,
                                                                       3691, 8667, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13672, 0, 3,
                                                                       13537, 8587, 13552, 3703,
                                                                       3721, 8677, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13717, 0, 3,
                                                                       13552, 8597, 13567, 3721,
                                                                       3739, 8707, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13762, 0, 3,
                                                                       13567, 8607, 13582, 3739,
                                                                       3757, 8737, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13807, 0, 3,
                                                                       13582, 8617, 13597, 3757,
                                                                       3775, 8767, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13852, 0, 3,
                                                                       13597, 8627, 13612, 3775,
                                                                       3793, 8797, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13897, 0, 3,
                                                                       13612, 8637, 13627, 3793,
                                                                       3811, 8827, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13942, 0, 3,
                                                                       13627, 8647, 13642, 3811,
                                                                       3829, 8857, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13987, 0, 3,
                                                                       13642, 8657, 13657, 3829,
                                                                       3847, 8887, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14032, 0, 3,
                                                                       13672, 8677, 13717, 3883,
                                                                       3919, 8917, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14122, 0, 3,
                                                                       13717, 8707, 13762, 3919,
                                                                       3955, 8977, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14212, 0, 3,
                                                                       13762, 8737, 13807, 3955,
                                                                       3991, 9037, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14302, 0, 3,
                                                                       13807, 8767, 13852, 3991,
                                                                       4027, 9097, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14392, 0, 3,
                                                                       13852, 8797, 13897, 4027,
                                                                       4063, 9157, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14482, 0, 3,
                                                                       13897, 8827, 13942, 4063,
                                                                       4099, 9217, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14572, 0, 3,
                                                                       13942, 8857, 13987, 4099,
                                                                       4135, 9277, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 14662, 0, 3,
                                                                       14032, 8917, 14122, 4207,
                                                                       4267, 9337, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 14812, 0, 3,
                                                                       14122, 8977, 14212, 4267,
                                                                       4327, 9437, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 14962, 0, 3,
                                                                       14212, 9037, 14302, 4327,
                                                                       4387, 9537, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 15112, 0, 3,
                                                                       14302, 9097, 14392, 4387,
                                                                       4447, 9637, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 15262, 0, 3,
                                                                       14392, 9157, 14482, 4447,
                                                                       4507, 9737, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 15412, 0, 3,
                                                                       14482, 9217, 14572, 4507,
                                                                       4567, 9837, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 15562, 0, 3,
                                                                       14662, 9337, 14812, 4687,
                                                                       4777, 9937, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 15787, 0, 3,
                                                                       14812, 9437, 14962, 4777,
                                                                       4867, 10087, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 16012, 0, 3,
                                                                       14962, 9537, 15112, 4867,
                                                                       4957, 10237, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 16237, 0, 3,
                                                                       15112, 9637, 15262, 4957,
                                                                       5047, 10387, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 16462, 0, 3,
                                                                       15262, 9737, 15412, 5047,
                                                                       5137, 10537, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 16687, 0, 3,
                                                                       15562, 9937, 15787, 5317,
                                                                       5443, 10687, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 17002, 0, 3,
                                                                       15787, 10087, 16012, 5443,
                                                                       5569, 10897, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 17317, 0, 3,
                                                                       16012, 10237, 16237, 5569,
                                                                       5695, 11107, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 17632, 0, 3,
                                                                       16237, 10387, 16462, 5695,
                                                                       5821, 11317, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 17947, 0, 3,
                                                                       16687, 10687, 17002, 6073,
                                                                       6241, 11527, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 18367, 0, 3,
                                                                       17002, 10897, 17317, 6241,
                                                                       6409, 11807, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 18787, 0, 3,
                                                                       17317, 11107, 17632, 6409,
                                                                       6577, 12087, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 19207, 0, 3,
                                                                       17947, 11527, 18367, 6913,
                                                                       7129, 12367, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 19747, 0, 3,
                                                                       18367, 11807, 18787, 7129,
                                                                       7345, 12727, ncols, gamma,
                                                                       p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 20287, 0, 3,
                                                                       19207, 12367, 19747, 7777,
                                                                       8047, 13087, ncols, gamma,
                                                                       p, q);

                    simdgeo::geom_f_x(buffer, 20962, 14032, 15562, 1, 15, ncols, beta);

                    simdgeo::geom_f_y(buffer, 21112, 14032, 15562, 1, 15, ncols, beta);

                    simdgeo::geom_f_z(buffer, 21262, 14032, 15562, 1, 15, ncols, beta);

                    simdgeo::geom_g_x(buffer, 21412, 14662, 16687, 1, 15, ncols, beta);

                    simdgeo::geom_g_y(buffer, 21637, 14662, 16687, 1, 15, ncols, beta);

                    simdgeo::geom_g_z(buffer, 21862, 14662, 16687, 1, 15, ncols, beta);

                    simdgeo::geom_h_x(buffer, 22087, 15562, 17947, 1, 15, ncols, beta);

                    simdgeo::geom_h_y(buffer, 22402, 15562, 17947, 1, 15, ncols, beta);

                    simdgeo::geom_h_z(buffer, 22717, 15562, 17947, 1, 15, ncols, beta);

                    simdgeo::geom_i_x(buffer, 23032, 16687, 19207, 1, 15, ncols, beta);

                    simdgeo::geom_i_y(buffer, 23452, 16687, 19207, 1, 15, ncols, beta);

                    simdgeo::geom_i_z(buffer, 23872, 16687, 19207, 1, 15, ncols, beta);

                    simdgeo::geom_k_x(buffer, 24292, 17947, 20287, 1, 15, ncols, beta);

                    simdgeo::geom_k_y(buffer, 24832, 17947, 20287, 1, 15, ncols, beta);

                    simdgeo::geom_k_z(buffer, 25372, 17947, 20287, 1, 15, ncols, beta);

                    simdfunc::contract_primitives(buffer, 25912, 20962, 150, ncols);

                    simdfunc::contract_primitives(buffer, 26152, 21112, 150, ncols);

                    simdfunc::contract_primitives(buffer, 26392, 21262, 150, ncols);

                    simdfunc::contract_primitives(buffer, 26632, 14662, 150, ncols);

                    simdfunc::contract_primitives(buffer, 26872, 21412, 225, ncols);

                    simdfunc::contract_primitives(buffer, 27232, 21637, 225, ncols);

                    simdfunc::contract_primitives(buffer, 27592, 21862, 225, ncols);

                    simdfunc::contract_primitives(buffer, 27952, 15562, 225, ncols);

                    simdfunc::contract_primitives(buffer, 28312, 22087, 315, ncols);

                    simdfunc::contract_primitives(buffer, 28816, 22402, 315, ncols);

                    simdfunc::contract_primitives(buffer, 29320, 22717, 315, ncols);

                    simdfunc::contract_primitives(buffer, 29824, 16687, 315, ncols);

                    simdfunc::contract_primitives(buffer, 30328, 23032, 420, ncols);

                    simdfunc::contract_primitives(buffer, 31000, 23452, 420, ncols);

                    simdfunc::contract_primitives(buffer, 31672, 23872, 420, ncols);

                    simdfunc::contract_primitives(buffer, 32344, 17947, 420, ncols);

                    simdfunc::contract_primitives(buffer, 33016, 24292, 540, ncols);

                    simdfunc::contract_primitives(buffer, 33880, 24832, 540, ncols);

                    simdfunc::contract_primitives(buffer, 34744, 25372, 540, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 26062, 25912, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 26302, 26152, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 26542, 26392, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 26782, 26632, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 27097, 26872, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 27457, 27232, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 27817, 27592, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 28177, 27952, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 28627, 28312, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 29131, 28816, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 29635, 29320, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 30139, 29824, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 30748, 30328, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 31420, 31000, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 32092, 31672, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 32764, 32344, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 33556, 33016, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 34420, 33880, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 35284, 34744, 36, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 35608, 26062, 26782, 27097, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 35878, 26302, 26782, 27457, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 36148, 26542, 26782, 27817, 9,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 36418, 26782, 28177, 9, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 36688, 27097, 28177, 28627, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 37093, 27457, 28177, 29131, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 37498, 27817, 28177, 29635, 9,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 37903, 28177, 30139, 9, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 38308, 28627, 30139, 30748, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 38875, 29131, 30139, 31420, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 39442, 29635, 30139, 32092, 9,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 40009, 30139, 32764, 9, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 40576, 30748, 32764, 33556, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 41332, 31420, 32764, 34420, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 42088, 32092, 32764, 35284, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 42844, 35608, 36418, 36688, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 43384, 35878, 36418, 37093, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 43924, 36148, 36418, 37498, 9,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 44464, 36418, 37903, 9, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 45004, 36688, 37903, 38308, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 45814, 37093, 37903, 38875, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 46624, 37498, 37903, 39442, 9,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 47434, 37903, 40009, 9, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 48244, 38308, 40009, 40576, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 49378, 38875, 40009, 41332, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 50512, 39442, 40009, 42088, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 51646, 42844, 44464, 45004, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 52546, 43384, 44464, 45814, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 53446, 43924, 44464, 46624, 9,
                                          nmax);

        simdtrf::compute_hrr_ff(buffer, coordinates, 54346, 44464, 47434, 9, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 55246, 45004, 47434, 48244, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 56596, 45814, 47434, 49378, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 57946, 46624, 47434, 50512, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gf_out_of_second(buffer, coordinates, 59296, 51646, 54346,
                                                        55246, 9, nmax);

        simdtrf::compute_hrr_geom_010y_gf_out_of_second(buffer, coordinates, 60646, 52546, 54346,
                                                        56596, 9, nmax);

        simdtrf::compute_hrr_geom_010z_gf_out_of_second(buffer, coordinates, 61996, 53446, 54346,
                                                        57946, 9, nmax);

        simdtrf::transform_f_inner(buffer, 63346, 59296, 15, 9, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 63346, 63, nmax);

        simdtrf::transform_f_inner(buffer, 63346, 60646, 15, 9, nmax);

        simdtrf::transform_g_outer(values + 567 * nvalues + n * npairs, nvalues, buffer, 63346,
                                   63, nmax);

        simdtrf::transform_f_inner(buffer, 63346, 61996, 15, 9, nmax);

        simdtrf::transform_g_outer(values + 1134 * nvalues + n * npairs, nvalues, buffer, 63346,
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
