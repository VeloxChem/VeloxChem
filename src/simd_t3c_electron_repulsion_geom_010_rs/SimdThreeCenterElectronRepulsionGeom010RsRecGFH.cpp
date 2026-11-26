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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGFH.hpp"

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
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_gfh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_gfh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 179927, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 4158 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 179927, 85112, 25460, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10, 11, 12, 13},
                                                            ncols, fj, i * nprim_b + j, fq,
                                                            omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 20, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13}, ncols,
                                                        fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 61, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 64, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 67, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 70, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 73, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 76, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 79, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 82, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 85, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 88, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 91, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 94, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 97, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 100, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 103, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 106, 0, 3, 7, 8,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 112, 0, 3, 8, 9,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 118, 0, 3, 9, 10,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 124, 0, 3, 10, 11,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 130, 0, 3, 11, 12,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 136, 0, 3, 12, 13,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 142, 0, 3, 13, 14,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 148, 0, 3, 14, 15,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 154, 0, 3, 15, 16,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 160, 0, 3, 16, 17,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 166, 0, 3, 17, 18,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 172, 0, 3, 21, 22,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 178, 0, 3, 22, 23,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 184, 0, 3, 23, 24,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 190, 0, 3, 24, 25,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 196, 0, 3, 25, 26,
                                                                       82, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 202, 0, 3, 26, 27,
                                                                       85, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 208, 0, 3, 27, 28,
                                                                       88, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 214, 0, 3, 28, 29,
                                                                       91, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 220, 0, 3, 29, 30,
                                                                       94, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 226, 0, 3, 30, 31,
                                                                       97, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 232, 0, 3, 31, 32,
                                                                       100, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 34, 37,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 248, 0, 3, 37, 40,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 40, 43,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 43, 46,
                                                                       124, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 46, 49,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 49, 52,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 52, 55,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 55, 58,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 58, 61,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 61, 64,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 70, 73,
                                                                       172, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 73, 76,
                                                                       178, 184, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 76, 79,
                                                                       184, 190, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 79, 82,
                                                                       190, 196, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 82, 85,
                                                                       196, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 85, 88,
                                                                       202, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 88, 91,
                                                                       208, 214, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 91, 94,
                                                                       214, 220, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 94, 97,
                                                                       220, 226, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 97,
                                                                       100, 226, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 106,
                                                                       112, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 453, 0, 3, 112,
                                                                       118, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 118,
                                                                       124, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 483, 0, 3, 124,
                                                                       130, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 130,
                                                                       136, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 513, 0, 3, 136,
                                                                       142, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 528, 0, 3, 142,
                                                                       148, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 543, 0, 3, 148,
                                                                       154, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 558, 0, 3, 154,
                                                                       160, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 573, 0, 3, 172,
                                                                       178, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 588, 0, 3, 178,
                                                                       184, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 603, 0, 3, 184,
                                                                       190, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 618, 0, 3, 190,
                                                                       196, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 633, 0, 3, 196,
                                                                       202, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 648, 0, 3, 202,
                                                                       208, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 663, 0, 3, 208,
                                                                       214, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 678, 0, 3, 214,
                                                                       220, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 693, 0, 3, 220,
                                                                       226, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 708, 0, 3, 238,
                                                                       248, 438, 453, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 729, 0, 3, 248,
                                                                       258, 453, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 750, 0, 3, 258,
                                                                       268, 468, 483, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 771, 0, 3, 268,
                                                                       278, 483, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 792, 0, 3, 278,
                                                                       288, 498, 513, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 813, 0, 3, 288,
                                                                       298, 513, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 834, 0, 3, 298,
                                                                       308, 528, 543, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 855, 0, 3, 308,
                                                                       318, 543, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 876, 0, 3, 338,
                                                                       348, 573, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 897, 0, 3, 348,
                                                                       358, 588, 603, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 918, 0, 3, 358,
                                                                       368, 603, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 939, 0, 3, 368,
                                                                       378, 618, 633, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 960, 0, 3, 378,
                                                                       388, 633, 648, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 981, 0, 3, 388,
                                                                       398, 648, 663, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 398,
                                                                       408, 663, 678, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1023, 0, 3, 408,
                                                                       418, 678, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 438,
                                                                       453, 708, 729, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 453,
                                                                       468, 729, 750, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 468,
                                                                       483, 750, 771, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 483,
                                                                       498, 771, 792, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 498,
                                                                       513, 792, 813, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 513,
                                                                       528, 813, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 528,
                                                                       543, 834, 855, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 573,
                                                                       588, 876, 897, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 588,
                                                                       603, 897, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 603,
                                                                       618, 918, 939, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 618,
                                                                       633, 939, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 633,
                                                                       648, 960, 981, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 648,
                                                                       663, 981, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 663,
                                                                       678, 1002, 1023, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 708,
                                                                       729, 1044, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1472, 0, 3, 729,
                                                                       750, 1072, 1100, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1508, 0, 3, 750,
                                                                       771, 1100, 1128, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1544, 0, 3, 771,
                                                                       792, 1128, 1156, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1580, 0, 3, 792,
                                                                       813, 1156, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1616, 0, 3, 813,
                                                                       834, 1184, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1652, 0, 3, 876,
                                                                       897, 1240, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 897,
                                                                       918, 1268, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1724, 0, 3, 918,
                                                                       939, 1296, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1760, 0, 3, 939,
                                                                       960, 1324, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1796, 0, 3, 960,
                                                                       981, 1352, 1380, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1832, 0, 3, 981,
                                                                       1002, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1868, 0, 3, 1044,
                                                                       1072, 1436, 1472, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1913, 0, 3, 1072,
                                                                       1100, 1472, 1508, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1958, 0, 3, 1100,
                                                                       1128, 1508, 1544, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2003, 0, 3, 1128,
                                                                       1156, 1544, 1580, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2048, 0, 3, 1156,
                                                                       1184, 1580, 1616, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2093, 0, 3, 1240,
                                                                       1268, 1652, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2138, 0, 3, 1268,
                                                                       1296, 1688, 1724, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2183, 0, 3, 1296,
                                                                       1324, 1724, 1760, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2228, 0, 3, 1324,
                                                                       1352, 1760, 1796, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2273, 0, 3, 1352,
                                                                       1380, 1796, 1832, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2318, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2321, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2324, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2327, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2330, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2333, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2336, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2339, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2342, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2345, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2348, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2351, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2354, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2357, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2360, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2363, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2366, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2369, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2372, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2375, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2378, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2381, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2384, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2387, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2390, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2393, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2396, 3, 9, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2405, 3, 10, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2414, 3, 11, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2423, 3, 12, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2432, 3, 13, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2441, 3, 14, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2450, 3, 15, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2459, 3, 16, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2468, 3, 17, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2477, 3, 18, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2486, 3, 23, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2495, 3, 24, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2504, 3, 25, 82,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2513, 3, 26, 85,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2522, 3, 27, 88,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2531, 3, 28, 91,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2540, 3, 29, 94,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2549, 3, 30, 97,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2558, 3, 31, 100,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2567, 3, 32, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2576, 3, 34, 106,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2594, 3, 37, 112,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2612, 3, 40, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2630, 3, 43, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2648, 3, 46, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2666, 3, 49, 136,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2684, 3, 52, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2702, 3, 55, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2720, 3, 58, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2738, 3, 61, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2756, 3, 64, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2774, 3, 70, 172,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2792, 3, 73, 178,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2810, 3, 76, 184,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2828, 3, 79, 190,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2846, 3, 82, 196,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2864, 3, 85, 202,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2882, 3, 88, 208,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2900, 3, 91, 214,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2918, 3, 94, 220,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2936, 3, 97, 226,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2954, 3, 100, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2972, 3, 106, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3002, 3, 112, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3032, 3, 118, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3062, 3, 124, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3092, 3, 130, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3122, 3, 136, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3152, 3, 142, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3182, 3, 148, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3212, 3, 154, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3242, 3, 160, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3272, 3, 172, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3302, 3, 178, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3332, 3, 184, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3362, 3, 190, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3392, 3, 196, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3422, 3, 202, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3452, 3, 208, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3482, 3, 214, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3512, 3, 220, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3542, 3, 226, 428,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3572, 3, 238, 438,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3617, 3, 248, 453,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3662, 3, 258, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3707, 3, 268, 483,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3752, 3, 278, 498,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3797, 3, 288, 513,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3842, 3, 298, 528,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3887, 3, 308, 543,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3932, 3, 318, 558,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3977, 3, 338, 573,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4022, 3, 348, 588,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4067, 3, 358, 603,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4112, 3, 368, 618,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4157, 3, 378, 633,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4202, 3, 388, 648,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4247, 3, 398, 663,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4292, 3, 408, 678,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4337, 3, 418, 693,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4382, 3, 438, 708,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4445, 3, 453, 729,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4508, 3, 468, 750,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4571, 3, 483, 771,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4634, 3, 498, 792,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4697, 3, 513, 813,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4760, 3, 528, 834,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4823, 3, 543, 855,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4886, 3, 573, 876,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4949, 3, 588, 897,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5012, 3, 603, 918,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5075, 3, 618, 939,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5138, 3, 633, 960,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5201, 3, 648, 981,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5264, 3, 663,
                                                                       1002, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5327, 3, 678,
                                                                       1023, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5390, 3, 708,
                                                                       1044, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5474, 3, 729,
                                                                       1072, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5558, 3, 750,
                                                                       1100, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5642, 3, 771,
                                                                       1128, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5726, 3, 792,
                                                                       1156, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5810, 3, 813,
                                                                       1184, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5894, 3, 834,
                                                                       1212, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5978, 3, 876,
                                                                       1240, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6062, 3, 897,
                                                                       1268, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6146, 3, 918,
                                                                       1296, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6230, 3, 939,
                                                                       1324, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6314, 3, 960,
                                                                       1352, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6398, 3, 981,
                                                                       1380, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6482, 3, 1002,
                                                                       1408, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6566, 3, 1044,
                                                                       1436, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6674, 3, 1072,
                                                                       1472, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6782, 3, 1100,
                                                                       1508, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6890, 3, 1128,
                                                                       1544, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6998, 3, 1156,
                                                                       1580, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7106, 3, 1184,
                                                                       1616, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7214, 3, 1240,
                                                                       1652, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7322, 3, 1268,
                                                                       1688, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7430, 3, 1296,
                                                                       1724, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7538, 3, 1324,
                                                                       1760, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7646, 3, 1352,
                                                                       1796, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7754, 3, 1380,
                                                                       1832, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7862, 3, 1436,
                                                                       1868, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7997, 3, 1472,
                                                                       1913, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8132, 3, 1508,
                                                                       1958, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8267, 3, 1544,
                                                                       2003, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8402, 3, 1580,
                                                                       2048, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8537, 3, 1652,
                                                                       2093, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8672, 3, 1688,
                                                                       2138, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8807, 3, 1724,
                                                                       2183, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8942, 3, 1760,
                                                                       2228, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9077, 3, 1796,
                                                                       2273, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9212, 3, 7, 8,
                                                                       2324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9218, 3, 8, 9,
                                                                       2327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9224, 3, 9, 10,
                                                                       2330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9230, 3, 10, 11,
                                                                       2333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9236, 3, 11, 12,
                                                                       2336, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9242, 3, 12, 13,
                                                                       2339, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9248, 3, 13, 14,
                                                                       2342, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9254, 3, 14, 15,
                                                                       2345, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9260, 3, 15, 16,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9266, 3, 16, 17,
                                                                       2351, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9272, 3, 17, 18,
                                                                       2354, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9278, 3, 21, 22,
                                                                       2363, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9284, 3, 22, 23,
                                                                       2366, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9290, 3, 23, 24,
                                                                       2369, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9296, 3, 24, 25,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9302, 3, 25, 26,
                                                                       2375, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9308, 3, 26, 27,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9314, 3, 27, 28,
                                                                       2381, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9320, 3, 28, 29,
                                                                       2384, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9326, 3, 29, 30,
                                                                       2387, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9332, 3, 30, 31,
                                                                       2390, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9338, 3, 31, 32,
                                                                       2393, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9344, 0, 3, 9212,
                                                                       2324, 9218, 2396, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9362, 0, 3, 9218,
                                                                       2327, 9224, 2405, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9380, 0, 3, 9224,
                                                                       2330, 9230, 2414, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9398, 0, 3, 9230,
                                                                       2333, 9236, 2423, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9416, 0, 3, 9236,
                                                                       2336, 9242, 2432, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9434, 0, 3, 9242,
                                                                       2339, 9248, 2441, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9452, 0, 3, 9248,
                                                                       2342, 9254, 2450, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9470, 0, 3, 9254,
                                                                       2345, 9260, 2459, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9488, 0, 3, 9260,
                                                                       2348, 9266, 2468, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9506, 0, 3, 9266,
                                                                       2351, 9272, 2477, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9524, 0, 3, 9278,
                                                                       2363, 9284, 2486, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9542, 0, 3, 9284,
                                                                       2366, 9290, 2495, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9560, 0, 3, 9290,
                                                                       2369, 9296, 2504, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9578, 0, 3, 9296,
                                                                       2372, 9302, 2513, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9596, 0, 3, 9302,
                                                                       2375, 9308, 2522, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9614, 0, 3, 9308,
                                                                       2378, 9314, 2531, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9632, 0, 3, 9314,
                                                                       2381, 9320, 2540, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9650, 0, 3, 9320,
                                                                       2384, 9326, 2549, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9668, 0, 3, 9326,
                                                                       2387, 9332, 2558, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9686, 0, 3, 9332,
                                                                       2390, 9338, 2567, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9704, 0, 3, 9344,
                                                                       2396, 9362, 106, 112,
                                                                       2612, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9740, 0, 3, 9362,
                                                                       2405, 9380, 112, 118,
                                                                       2630, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9776, 0, 3, 9380,
                                                                       2414, 9398, 118, 124,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9812, 0, 3, 9398,
                                                                       2423, 9416, 124, 130,
                                                                       2666, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9848, 0, 3, 9416,
                                                                       2432, 9434, 130, 136,
                                                                       2684, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9884, 0, 3, 9434,
                                                                       2441, 9452, 136, 142,
                                                                       2702, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9920, 0, 3, 9452,
                                                                       2450, 9470, 142, 148,
                                                                       2720, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9956, 0, 3, 9470,
                                                                       2459, 9488, 148, 154,
                                                                       2738, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9992, 0, 3, 9488,
                                                                       2468, 9506, 154, 160,
                                                                       2756, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10028, 0, 3, 9524,
                                                                       2486, 9542, 172, 178,
                                                                       2810, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10064, 0, 3, 9542,
                                                                       2495, 9560, 178, 184,
                                                                       2828, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10100, 0, 3, 9560,
                                                                       2504, 9578, 184, 190,
                                                                       2846, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10136, 0, 3, 9578,
                                                                       2513, 9596, 190, 196,
                                                                       2864, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10172, 0, 3, 9596,
                                                                       2522, 9614, 196, 202,
                                                                       2882, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10208, 0, 3, 9614,
                                                                       2531, 9632, 202, 208,
                                                                       2900, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10244, 0, 3, 9632,
                                                                       2540, 9650, 208, 214,
                                                                       2918, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10280, 0, 3, 9650,
                                                                       2549, 9668, 214, 220,
                                                                       2936, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10316, 0, 3, 9668,
                                                                       2558, 9686, 220, 226,
                                                                       2954, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10352, 0, 3, 9704,
                                                                       2612, 9740, 238, 248,
                                                                       3032, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10412, 0, 3, 9740,
                                                                       2630, 9776, 248, 258,
                                                                       3062, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10472, 0, 3, 9776,
                                                                       2648, 9812, 258, 268,
                                                                       3092, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10532, 0, 3, 9812,
                                                                       2666, 9848, 268, 278,
                                                                       3122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10592, 0, 3, 9848,
                                                                       2684, 9884, 278, 288,
                                                                       3152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10652, 0, 3, 9884,
                                                                       2702, 9920, 288, 298,
                                                                       3182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10712, 0, 3, 9920,
                                                                       2720, 9956, 298, 308,
                                                                       3212, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10772, 0, 3, 9956,
                                                                       2738, 9992, 308, 318,
                                                                       3242, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10832, 0, 3,
                                                                       10028, 2810, 10064, 338,
                                                                       348, 3332, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10892, 0, 3,
                                                                       10064, 2828, 10100, 348,
                                                                       358, 3362, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10952, 0, 3,
                                                                       10100, 2846, 10136, 358,
                                                                       368, 3392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11012, 0, 3,
                                                                       10136, 2864, 10172, 368,
                                                                       378, 3422, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11072, 0, 3,
                                                                       10172, 2882, 10208, 378,
                                                                       388, 3452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11132, 0, 3,
                                                                       10208, 2900, 10244, 388,
                                                                       398, 3482, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11192, 0, 3,
                                                                       10244, 2918, 10280, 398,
                                                                       408, 3512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11252, 0, 3,
                                                                       10280, 2936, 10316, 408,
                                                                       418, 3542, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11312, 0, 3,
                                                                       10352, 3032, 10412, 438,
                                                                       453, 3662, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11402, 0, 3,
                                                                       10412, 3062, 10472, 453,
                                                                       468, 3707, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11492, 0, 3,
                                                                       10472, 3092, 10532, 468,
                                                                       483, 3752, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11582, 0, 3,
                                                                       10532, 3122, 10592, 483,
                                                                       498, 3797, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11672, 0, 3,
                                                                       10592, 3152, 10652, 498,
                                                                       513, 3842, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11762, 0, 3,
                                                                       10652, 3182, 10712, 513,
                                                                       528, 3887, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11852, 0, 3,
                                                                       10712, 3212, 10772, 528,
                                                                       543, 3932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11942, 0, 3,
                                                                       10832, 3332, 10892, 573,
                                                                       588, 4067, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12032, 0, 3,
                                                                       10892, 3362, 10952, 588,
                                                                       603, 4112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12122, 0, 3,
                                                                       10952, 3392, 11012, 603,
                                                                       618, 4157, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12212, 0, 3,
                                                                       11012, 3422, 11072, 618,
                                                                       633, 4202, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12302, 0, 3,
                                                                       11072, 3452, 11132, 633,
                                                                       648, 4247, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12392, 0, 3,
                                                                       11132, 3482, 11192, 648,
                                                                       663, 4292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12482, 0, 3,
                                                                       11192, 3512, 11252, 663,
                                                                       678, 4337, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12572, 0, 3,
                                                                       11312, 3662, 11402, 708,
                                                                       729, 4508, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12698, 0, 3,
                                                                       11402, 3707, 11492, 729,
                                                                       750, 4571, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12824, 0, 3,
                                                                       11492, 3752, 11582, 750,
                                                                       771, 4634, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12950, 0, 3,
                                                                       11582, 3797, 11672, 771,
                                                                       792, 4697, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13076, 0, 3,
                                                                       11672, 3842, 11762, 792,
                                                                       813, 4760, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13202, 0, 3,
                                                                       11762, 3887, 11852, 813,
                                                                       834, 4823, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13328, 0, 3,
                                                                       11942, 4067, 12032, 876,
                                                                       897, 5012, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13454, 0, 3,
                                                                       12032, 4112, 12122, 897,
                                                                       918, 5075, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13580, 0, 3,
                                                                       12122, 4157, 12212, 918,
                                                                       939, 5138, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13706, 0, 3,
                                                                       12212, 4202, 12302, 939,
                                                                       960, 5201, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13832, 0, 3,
                                                                       12302, 4247, 12392, 960,
                                                                       981, 5264, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13958, 0, 3,
                                                                       12392, 4292, 12482, 981,
                                                                       1002, 5327, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14084, 0, 3,
                                                                       12572, 4508, 12698, 1044,
                                                                       1072, 5558, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14252, 0, 3,
                                                                       12698, 4571, 12824, 1072,
                                                                       1100, 5642, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14420, 0, 3,
                                                                       12824, 4634, 12950, 1100,
                                                                       1128, 5726, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14588, 0, 3,
                                                                       12950, 4697, 13076, 1128,
                                                                       1156, 5810, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14756, 0, 3,
                                                                       13076, 4760, 13202, 1156,
                                                                       1184, 5894, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14924, 0, 3,
                                                                       13328, 5012, 13454, 1240,
                                                                       1268, 6146, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15092, 0, 3,
                                                                       13454, 5075, 13580, 1268,
                                                                       1296, 6230, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15260, 0, 3,
                                                                       13580, 5138, 13706, 1296,
                                                                       1324, 6314, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15428, 0, 3,
                                                                       13706, 5201, 13832, 1324,
                                                                       1352, 6398, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15596, 0, 3,
                                                                       13832, 5264, 13958, 1352,
                                                                       1380, 6482, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15764, 0, 3,
                                                                       14084, 5558, 14252, 1436,
                                                                       1472, 6782, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15980, 0, 3,
                                                                       14252, 5642, 14420, 1472,
                                                                       1508, 6890, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16196, 0, 3,
                                                                       14420, 5726, 14588, 1508,
                                                                       1544, 6998, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16412, 0, 3,
                                                                       14588, 5810, 14756, 1544,
                                                                       1580, 7106, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16628, 0, 3,
                                                                       14924, 6146, 15092, 1652,
                                                                       1688, 7430, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16844, 0, 3,
                                                                       15092, 6230, 15260, 1688,
                                                                       1724, 7538, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17060, 0, 3,
                                                                       15260, 6314, 15428, 1724,
                                                                       1760, 7646, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17276, 0, 3,
                                                                       15428, 6398, 15596, 1760,
                                                                       1796, 7754, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17492, 0, 3,
                                                                       15764, 6782, 15980, 1868,
                                                                       1913, 8132, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17762, 0, 3,
                                                                       15980, 6890, 16196, 1913,
                                                                       1958, 8267, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18032, 0, 3,
                                                                       16196, 6998, 16412, 1958,
                                                                       2003, 8402, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18302, 0, 3,
                                                                       16628, 7430, 16844, 2093,
                                                                       2138, 8807, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18572, 0, 3,
                                                                       16844, 7538, 17060, 2138,
                                                                       2183, 8942, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18842, 0, 3,
                                                                       17060, 7646, 17276, 2183,
                                                                       2228, 9077, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19112, 3, 2318,
                                                                       2321, 9212, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19122, 3, 2321,
                                                                       2324, 9218, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19132, 3, 2324,
                                                                       2327, 9224, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19142, 3, 2327,
                                                                       2330, 9230, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19152, 3, 2330,
                                                                       2333, 9236, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19162, 3, 2333,
                                                                       2336, 9242, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19172, 3, 2336,
                                                                       2339, 9248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19182, 3, 2339,
                                                                       2342, 9254, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19192, 3, 2342,
                                                                       2345, 9260, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19202, 3, 2345,
                                                                       2348, 9266, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19212, 3, 2348,
                                                                       2351, 9272, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19222, 3, 2357,
                                                                       2360, 9278, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19232, 3, 2360,
                                                                       2363, 9284, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19242, 3, 2363,
                                                                       2366, 9290, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19252, 3, 2366,
                                                                       2369, 9296, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19262, 3, 2369,
                                                                       2372, 9302, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19272, 3, 2372,
                                                                       2375, 9308, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19282, 3, 2375,
                                                                       2378, 9314, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19292, 3, 2378,
                                                                       2381, 9320, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19302, 3, 2381,
                                                                       2384, 9326, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19312, 3, 2384,
                                                                       2387, 9332, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19322, 3, 2387,
                                                                       2390, 9338, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19332, 0, 3,
                                                                       19112, 9212, 19122, 9344,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19362, 0, 3,
                                                                       19122, 9218, 19132, 9362,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19392, 0, 3,
                                                                       19132, 9224, 19142, 9380,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19422, 0, 3,
                                                                       19142, 9230, 19152, 9398,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19452, 0, 3,
                                                                       19152, 9236, 19162, 9416,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19482, 0, 3,
                                                                       19162, 9242, 19172, 9434,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19512, 0, 3,
                                                                       19172, 9248, 19182, 9452,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19542, 0, 3,
                                                                       19182, 9254, 19192, 9470,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19572, 0, 3,
                                                                       19192, 9260, 19202, 9488,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19602, 0, 3,
                                                                       19202, 9266, 19212, 9506,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19632, 0, 3,
                                                                       19222, 9278, 19232, 9524,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19662, 0, 3,
                                                                       19232, 9284, 19242, 9542,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19692, 0, 3,
                                                                       19242, 9290, 19252, 9560,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19722, 0, 3,
                                                                       19252, 9296, 19262, 9578,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19752, 0, 3,
                                                                       19262, 9302, 19272, 9596,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19782, 0, 3,
                                                                       19272, 9308, 19282, 9614,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19812, 0, 3,
                                                                       19282, 9314, 19292, 9632,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19842, 0, 3,
                                                                       19292, 9320, 19302, 9650,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19872, 0, 3,
                                                                       19302, 9326, 19312, 9668,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 19902, 0, 3,
                                                                       19312, 9332, 19322, 9686,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19932, 0, 3,
                                                                       19332, 9344, 19362, 2576,
                                                                       2594, 9704, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19992, 0, 3,
                                                                       19362, 9362, 19392, 2594,
                                                                       2612, 9740, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20052, 0, 3,
                                                                       19392, 9380, 19422, 2612,
                                                                       2630, 9776, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20112, 0, 3,
                                                                       19422, 9398, 19452, 2630,
                                                                       2648, 9812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20172, 0, 3,
                                                                       19452, 9416, 19482, 2648,
                                                                       2666, 9848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20232, 0, 3,
                                                                       19482, 9434, 19512, 2666,
                                                                       2684, 9884, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20292, 0, 3,
                                                                       19512, 9452, 19542, 2684,
                                                                       2702, 9920, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20352, 0, 3,
                                                                       19542, 9470, 19572, 2702,
                                                                       2720, 9956, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20412, 0, 3,
                                                                       19572, 9488, 19602, 2720,
                                                                       2738, 9992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20472, 0, 3,
                                                                       19632, 9524, 19662, 2774,
                                                                       2792, 10028, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20532, 0, 3,
                                                                       19662, 9542, 19692, 2792,
                                                                       2810, 10064, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20592, 0, 3,
                                                                       19692, 9560, 19722, 2810,
                                                                       2828, 10100, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20652, 0, 3,
                                                                       19722, 9578, 19752, 2828,
                                                                       2846, 10136, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20712, 0, 3,
                                                                       19752, 9596, 19782, 2846,
                                                                       2864, 10172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20772, 0, 3,
                                                                       19782, 9614, 19812, 2864,
                                                                       2882, 10208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20832, 0, 3,
                                                                       19812, 9632, 19842, 2882,
                                                                       2900, 10244, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20892, 0, 3,
                                                                       19842, 9650, 19872, 2900,
                                                                       2918, 10280, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20952, 0, 3,
                                                                       19872, 9668, 19902, 2918,
                                                                       2936, 10316, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21012, 0, 3,
                                                                       19932, 9704, 19992, 2972,
                                                                       3002, 10352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21112, 0, 3,
                                                                       19992, 9740, 20052, 3002,
                                                                       3032, 10412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21212, 0, 3,
                                                                       20052, 9776, 20112, 3032,
                                                                       3062, 10472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21312, 0, 3,
                                                                       20112, 9812, 20172, 3062,
                                                                       3092, 10532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21412, 0, 3,
                                                                       20172, 9848, 20232, 3092,
                                                                       3122, 10592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21512, 0, 3,
                                                                       20232, 9884, 20292, 3122,
                                                                       3152, 10652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21612, 0, 3,
                                                                       20292, 9920, 20352, 3152,
                                                                       3182, 10712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21712, 0, 3,
                                                                       20352, 9956, 20412, 3182,
                                                                       3212, 10772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21812, 0, 3,
                                                                       20472, 10028, 20532, 3272,
                                                                       3302, 10832, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21912, 0, 3,
                                                                       20532, 10064, 20592, 3302,
                                                                       3332, 10892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22012, 0, 3,
                                                                       20592, 10100, 20652, 3332,
                                                                       3362, 10952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22112, 0, 3,
                                                                       20652, 10136, 20712, 3362,
                                                                       3392, 11012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22212, 0, 3,
                                                                       20712, 10172, 20772, 3392,
                                                                       3422, 11072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22312, 0, 3,
                                                                       20772, 10208, 20832, 3422,
                                                                       3452, 11132, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22412, 0, 3,
                                                                       20832, 10244, 20892, 3452,
                                                                       3482, 11192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22512, 0, 3,
                                                                       20892, 10280, 20952, 3482,
                                                                       3512, 11252, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 22612, 0, 3,
                                                                       21012, 10352, 21112, 3572,
                                                                       3617, 11312, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 22762, 0, 3,
                                                                       21112, 10412, 21212, 3617,
                                                                       3662, 11402, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 22912, 0, 3,
                                                                       21212, 10472, 21312, 3662,
                                                                       3707, 11492, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23062, 0, 3,
                                                                       21312, 10532, 21412, 3707,
                                                                       3752, 11582, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23212, 0, 3,
                                                                       21412, 10592, 21512, 3752,
                                                                       3797, 11672, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23362, 0, 3,
                                                                       21512, 10652, 21612, 3797,
                                                                       3842, 11762, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23512, 0, 3,
                                                                       21612, 10712, 21712, 3842,
                                                                       3887, 11852, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23662, 0, 3,
                                                                       21812, 10832, 21912, 3977,
                                                                       4022, 11942, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23812, 0, 3,
                                                                       21912, 10892, 22012, 4022,
                                                                       4067, 12032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23962, 0, 3,
                                                                       22012, 10952, 22112, 4067,
                                                                       4112, 12122, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24112, 0, 3,
                                                                       22112, 11012, 22212, 4112,
                                                                       4157, 12212, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24262, 0, 3,
                                                                       22212, 11072, 22312, 4157,
                                                                       4202, 12302, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24412, 0, 3,
                                                                       22312, 11132, 22412, 4202,
                                                                       4247, 12392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24562, 0, 3,
                                                                       22412, 11192, 22512, 4247,
                                                                       4292, 12482, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 24712, 0, 3,
                                                                       22612, 11312, 22762, 4382,
                                                                       4445, 12572, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 24922, 0, 3,
                                                                       22762, 11402, 22912, 4445,
                                                                       4508, 12698, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25132, 0, 3,
                                                                       22912, 11492, 23062, 4508,
                                                                       4571, 12824, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25342, 0, 3,
                                                                       23062, 11582, 23212, 4571,
                                                                       4634, 12950, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25552, 0, 3,
                                                                       23212, 11672, 23362, 4634,
                                                                       4697, 13076, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25762, 0, 3,
                                                                       23362, 11762, 23512, 4697,
                                                                       4760, 13202, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25972, 0, 3,
                                                                       23662, 11942, 23812, 4886,
                                                                       4949, 13328, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26182, 0, 3,
                                                                       23812, 12032, 23962, 4949,
                                                                       5012, 13454, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26392, 0, 3,
                                                                       23962, 12122, 24112, 5012,
                                                                       5075, 13580, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26602, 0, 3,
                                                                       24112, 12212, 24262, 5075,
                                                                       5138, 13706, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26812, 0, 3,
                                                                       24262, 12302, 24412, 5138,
                                                                       5201, 13832, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 27022, 0, 3,
                                                                       24412, 12392, 24562, 5201,
                                                                       5264, 13958, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27232, 0, 3,
                                                                       24712, 12572, 24922, 5390,
                                                                       5474, 14084, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27512, 0, 3,
                                                                       24922, 12698, 25132, 5474,
                                                                       5558, 14252, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27792, 0, 3,
                                                                       25132, 12824, 25342, 5558,
                                                                       5642, 14420, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28072, 0, 3,
                                                                       25342, 12950, 25552, 5642,
                                                                       5726, 14588, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28352, 0, 3,
                                                                       25552, 13076, 25762, 5726,
                                                                       5810, 14756, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28632, 0, 3,
                                                                       25972, 13328, 26182, 5978,
                                                                       6062, 14924, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28912, 0, 3,
                                                                       26182, 13454, 26392, 6062,
                                                                       6146, 15092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 29192, 0, 3,
                                                                       26392, 13580, 26602, 6146,
                                                                       6230, 15260, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 29472, 0, 3,
                                                                       26602, 13706, 26812, 6230,
                                                                       6314, 15428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 29752, 0, 3,
                                                                       26812, 13832, 27022, 6314,
                                                                       6398, 15596, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 30032, 0, 3,
                                                                       27232, 14084, 27512, 6566,
                                                                       6674, 15764, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 30392, 0, 3,
                                                                       27512, 14252, 27792, 6674,
                                                                       6782, 15980, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 30752, 0, 3,
                                                                       27792, 14420, 28072, 6782,
                                                                       6890, 16196, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 31112, 0, 3,
                                                                       28072, 14588, 28352, 6890,
                                                                       6998, 16412, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 31472, 0, 3,
                                                                       28632, 14924, 28912, 7214,
                                                                       7322, 16628, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 31832, 0, 3,
                                                                       28912, 15092, 29192, 7322,
                                                                       7430, 16844, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32192, 0, 3,
                                                                       29192, 15260, 29472, 7430,
                                                                       7538, 17060, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32552, 0, 3,
                                                                       29472, 15428, 29752, 7538,
                                                                       7646, 17276, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 32912, 0, 3,
                                                                       30032, 15764, 30392, 7862,
                                                                       7997, 17492, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 33362, 0, 3,
                                                                       30392, 15980, 30752, 7997,
                                                                       8132, 17762, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 33812, 0, 3,
                                                                       30752, 16196, 31112, 8132,
                                                                       8267, 18032, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 34262, 0, 3,
                                                                       31472, 16628, 31832, 8537,
                                                                       8672, 18302, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 34712, 0, 3,
                                                                       31832, 16844, 32192, 8672,
                                                                       8807, 18572, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 35162, 0, 3,
                                                                       32192, 17060, 32552, 8807,
                                                                       8942, 18842, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35612, 3, 9212,
                                                                       9218, 19132, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35627, 3, 9218,
                                                                       9224, 19142, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35642, 3, 9224,
                                                                       9230, 19152, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35657, 3, 9230,
                                                                       9236, 19162, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35672, 3, 9236,
                                                                       9242, 19172, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35687, 3, 9242,
                                                                       9248, 19182, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35702, 3, 9248,
                                                                       9254, 19192, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35717, 3, 9254,
                                                                       9260, 19202, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35732, 3, 9260,
                                                                       9266, 19212, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35747, 3, 9278,
                                                                       9284, 19242, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35762, 3, 9284,
                                                                       9290, 19252, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35777, 3, 9290,
                                                                       9296, 19262, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35792, 3, 9296,
                                                                       9302, 19272, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35807, 3, 9302,
                                                                       9308, 19282, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35822, 3, 9308,
                                                                       9314, 19292, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35837, 3, 9314,
                                                                       9320, 19302, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35852, 3, 9320,
                                                                       9326, 19312, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35867, 3, 9326,
                                                                       9332, 19322, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35882, 0, 3,
                                                                       35612, 19132, 35627, 9344,
                                                                       9362, 19392, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35927, 0, 3,
                                                                       35627, 19142, 35642, 9362,
                                                                       9380, 19422, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35972, 0, 3,
                                                                       35642, 19152, 35657, 9380,
                                                                       9398, 19452, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36017, 0, 3,
                                                                       35657, 19162, 35672, 9398,
                                                                       9416, 19482, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36062, 0, 3,
                                                                       35672, 19172, 35687, 9416,
                                                                       9434, 19512, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36107, 0, 3,
                                                                       35687, 19182, 35702, 9434,
                                                                       9452, 19542, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36152, 0, 3,
                                                                       35702, 19192, 35717, 9452,
                                                                       9470, 19572, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36197, 0, 3,
                                                                       35717, 19202, 35732, 9470,
                                                                       9488, 19602, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36242, 0, 3,
                                                                       35747, 19242, 35762, 9524,
                                                                       9542, 19692, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36287, 0, 3,
                                                                       35762, 19252, 35777, 9542,
                                                                       9560, 19722, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36332, 0, 3,
                                                                       35777, 19262, 35792, 9560,
                                                                       9578, 19752, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36377, 0, 3,
                                                                       35792, 19272, 35807, 9578,
                                                                       9596, 19782, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36422, 0, 3,
                                                                       35807, 19282, 35822, 9596,
                                                                       9614, 19812, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36467, 0, 3,
                                                                       35822, 19292, 35837, 9614,
                                                                       9632, 19842, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36512, 0, 3,
                                                                       35837, 19302, 35852, 9632,
                                                                       9650, 19872, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 36557, 0, 3,
                                                                       35852, 19312, 35867, 9650,
                                                                       9668, 19902, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36602, 0, 3,
                                                                       35882, 19392, 35927, 9704,
                                                                       9740, 20052, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36692, 0, 3,
                                                                       35927, 19422, 35972, 9740,
                                                                       9776, 20112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36782, 0, 3,
                                                                       35972, 19452, 36017, 9776,
                                                                       9812, 20172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36872, 0, 3,
                                                                       36017, 19482, 36062, 9812,
                                                                       9848, 20232, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36962, 0, 3,
                                                                       36062, 19512, 36107, 9848,
                                                                       9884, 20292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37052, 0, 3,
                                                                       36107, 19542, 36152, 9884,
                                                                       9920, 20352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37142, 0, 3,
                                                                       36152, 19572, 36197, 9920,
                                                                       9956, 20412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37232, 0, 3,
                                                                       36242, 19692, 36287,
                                                                       10028, 10064, 20592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37322, 0, 3,
                                                                       36287, 19722, 36332,
                                                                       10064, 10100, 20652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37412, 0, 3,
                                                                       36332, 19752, 36377,
                                                                       10100, 10136, 20712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37502, 0, 3,
                                                                       36377, 19782, 36422,
                                                                       10136, 10172, 20772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37592, 0, 3,
                                                                       36422, 19812, 36467,
                                                                       10172, 10208, 20832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37682, 0, 3,
                                                                       36467, 19842, 36512,
                                                                       10208, 10244, 20892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 37772, 0, 3,
                                                                       36512, 19872, 36557,
                                                                       10244, 10280, 20952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37862, 0, 3,
                                                                       36602, 20052, 36692,
                                                                       10352, 10412, 21212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38012, 0, 3,
                                                                       36692, 20112, 36782,
                                                                       10412, 10472, 21312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38162, 0, 3,
                                                                       36782, 20172, 36872,
                                                                       10472, 10532, 21412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38312, 0, 3,
                                                                       36872, 20232, 36962,
                                                                       10532, 10592, 21512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38462, 0, 3,
                                                                       36962, 20292, 37052,
                                                                       10592, 10652, 21612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38612, 0, 3,
                                                                       37052, 20352, 37142,
                                                                       10652, 10712, 21712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38762, 0, 3,
                                                                       37232, 20592, 37322,
                                                                       10832, 10892, 22012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38912, 0, 3,
                                                                       37322, 20652, 37412,
                                                                       10892, 10952, 22112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39062, 0, 3,
                                                                       37412, 20712, 37502,
                                                                       10952, 11012, 22212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39212, 0, 3,
                                                                       37502, 20772, 37592,
                                                                       11012, 11072, 22312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39362, 0, 3,
                                                                       37592, 20832, 37682,
                                                                       11072, 11132, 22412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39512, 0, 3,
                                                                       37682, 20892, 37772,
                                                                       11132, 11192, 22512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 39662, 0, 3,
                                                                       37862, 21212, 38012,
                                                                       11312, 11402, 22912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 39887, 0, 3,
                                                                       38012, 21312, 38162,
                                                                       11402, 11492, 23062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40112, 0, 3,
                                                                       38162, 21412, 38312,
                                                                       11492, 11582, 23212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40337, 0, 3,
                                                                       38312, 21512, 38462,
                                                                       11582, 11672, 23362,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40562, 0, 3,
                                                                       38462, 21612, 38612,
                                                                       11672, 11762, 23512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40787, 0, 3,
                                                                       38762, 22012, 38912,
                                                                       11942, 12032, 23962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41012, 0, 3,
                                                                       38912, 22112, 39062,
                                                                       12032, 12122, 24112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41237, 0, 3,
                                                                       39062, 22212, 39212,
                                                                       12122, 12212, 24262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41462, 0, 3,
                                                                       39212, 22312, 39362,
                                                                       12212, 12302, 24412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41687, 0, 3,
                                                                       39362, 22412, 39512,
                                                                       12302, 12392, 24562,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 41912, 0, 3,
                                                                       39662, 22912, 39887,
                                                                       12572, 12698, 25132,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 42227, 0, 3,
                                                                       39887, 23062, 40112,
                                                                       12698, 12824, 25342,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 42542, 0, 3,
                                                                       40112, 23212, 40337,
                                                                       12824, 12950, 25552,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 42857, 0, 3,
                                                                       40337, 23362, 40562,
                                                                       12950, 13076, 25762,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43172, 0, 3,
                                                                       40787, 23962, 41012,
                                                                       13328, 13454, 26392,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43487, 0, 3,
                                                                       41012, 24112, 41237,
                                                                       13454, 13580, 26602,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43802, 0, 3,
                                                                       41237, 24262, 41462,
                                                                       13580, 13706, 26812,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 44117, 0, 3,
                                                                       41462, 24412, 41687,
                                                                       13706, 13832, 27022,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 44432, 0, 3,
                                                                       41912, 25132, 42227,
                                                                       14084, 14252, 27792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 44852, 0, 3,
                                                                       42227, 25342, 42542,
                                                                       14252, 14420, 28072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 45272, 0, 3,
                                                                       42542, 25552, 42857,
                                                                       14420, 14588, 28352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 45692, 0, 3,
                                                                       43172, 26392, 43487,
                                                                       14924, 15092, 29192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 46112, 0, 3,
                                                                       43487, 26602, 43802,
                                                                       15092, 15260, 29472,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 46532, 0, 3,
                                                                       43802, 26812, 44117,
                                                                       15260, 15428, 29752,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 46952, 0, 3,
                                                                       44432, 27792, 44852,
                                                                       15764, 15980, 30752,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 47492, 0, 3,
                                                                       44852, 28072, 45272,
                                                                       15980, 16196, 31112,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 48032, 0, 3,
                                                                       45692, 29192, 46112,
                                                                       16628, 16844, 32192,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 48572, 0, 3,
                                                                       46112, 29472, 46532,
                                                                       16844, 17060, 32552,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 49112, 0, 3,
                                                                       46952, 30752, 47492,
                                                                       17492, 17762, 33812,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 49787, 0, 3,
                                                                       48032, 32192, 48572,
                                                                       18302, 18572, 35162,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50462, 3, 19112,
                                                                       19122, 35612, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50483, 3, 19122,
                                                                       19132, 35627, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50504, 3, 19132,
                                                                       19142, 35642, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50525, 3, 19142,
                                                                       19152, 35657, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50546, 3, 19152,
                                                                       19162, 35672, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50567, 3, 19162,
                                                                       19172, 35687, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50588, 3, 19172,
                                                                       19182, 35702, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50609, 3, 19182,
                                                                       19192, 35717, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50630, 3, 19192,
                                                                       19202, 35732, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50651, 3, 19222,
                                                                       19232, 35747, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50672, 3, 19232,
                                                                       19242, 35762, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50693, 3, 19242,
                                                                       19252, 35777, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50714, 3, 19252,
                                                                       19262, 35792, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50735, 3, 19262,
                                                                       19272, 35807, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50756, 3, 19272,
                                                                       19282, 35822, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50777, 3, 19282,
                                                                       19292, 35837, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50798, 3, 19292,
                                                                       19302, 35852, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50819, 3, 19302,
                                                                       19312, 35867, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 50840, 0, 3,
                                                                       50462, 35612, 50483,
                                                                       19332, 19362, 35882,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 50903, 0, 3,
                                                                       50483, 35627, 50504,
                                                                       19362, 19392, 35927,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 50966, 0, 3,
                                                                       50504, 35642, 50525,
                                                                       19392, 19422, 35972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51029, 0, 3,
                                                                       50525, 35657, 50546,
                                                                       19422, 19452, 36017,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51092, 0, 3,
                                                                       50546, 35672, 50567,
                                                                       19452, 19482, 36062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51155, 0, 3,
                                                                       50567, 35687, 50588,
                                                                       19482, 19512, 36107,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51218, 0, 3,
                                                                       50588, 35702, 50609,
                                                                       19512, 19542, 36152,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51281, 0, 3,
                                                                       50609, 35717, 50630,
                                                                       19542, 19572, 36197,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51344, 0, 3,
                                                                       50651, 35747, 50672,
                                                                       19632, 19662, 36242,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51407, 0, 3,
                                                                       50672, 35762, 50693,
                                                                       19662, 19692, 36287,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51470, 0, 3,
                                                                       50693, 35777, 50714,
                                                                       19692, 19722, 36332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51533, 0, 3,
                                                                       50714, 35792, 50735,
                                                                       19722, 19752, 36377,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51596, 0, 3,
                                                                       50735, 35807, 50756,
                                                                       19752, 19782, 36422,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51659, 0, 3,
                                                                       50756, 35822, 50777,
                                                                       19782, 19812, 36467,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51722, 0, 3,
                                                                       50777, 35837, 50798,
                                                                       19812, 19842, 36512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51785, 0, 3,
                                                                       50798, 35852, 50819,
                                                                       19842, 19872, 36557,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 51848, 0, 3,
                                                                       50840, 35882, 50903,
                                                                       19932, 19992, 36602,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 51974, 0, 3,
                                                                       50903, 35927, 50966,
                                                                       19992, 20052, 36692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52100, 0, 3,
                                                                       50966, 35972, 51029,
                                                                       20052, 20112, 36782,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52226, 0, 3,
                                                                       51029, 36017, 51092,
                                                                       20112, 20172, 36872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52352, 0, 3,
                                                                       51092, 36062, 51155,
                                                                       20172, 20232, 36962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52478, 0, 3,
                                                                       51155, 36107, 51218,
                                                                       20232, 20292, 37052,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52604, 0, 3,
                                                                       51218, 36152, 51281,
                                                                       20292, 20352, 37142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52730, 0, 3,
                                                                       51344, 36242, 51407,
                                                                       20472, 20532, 37232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52856, 0, 3,
                                                                       51407, 36287, 51470,
                                                                       20532, 20592, 37322,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52982, 0, 3,
                                                                       51470, 36332, 51533,
                                                                       20592, 20652, 37412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53108, 0, 3,
                                                                       51533, 36377, 51596,
                                                                       20652, 20712, 37502,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53234, 0, 3,
                                                                       51596, 36422, 51659,
                                                                       20712, 20772, 37592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53360, 0, 3,
                                                                       51659, 36467, 51722,
                                                                       20772, 20832, 37682,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53486, 0, 3,
                                                                       51722, 36512, 51785,
                                                                       20832, 20892, 37772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 53612, 0, 3,
                                                                       51848, 36602, 51974,
                                                                       21012, 21112, 37862,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 53822, 0, 3,
                                                                       51974, 36692, 52100,
                                                                       21112, 21212, 38012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54032, 0, 3,
                                                                       52100, 36782, 52226,
                                                                       21212, 21312, 38162,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54242, 0, 3,
                                                                       52226, 36872, 52352,
                                                                       21312, 21412, 38312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54452, 0, 3,
                                                                       52352, 36962, 52478,
                                                                       21412, 21512, 38462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54662, 0, 3,
                                                                       52478, 37052, 52604,
                                                                       21512, 21612, 38612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54872, 0, 3,
                                                                       52730, 37232, 52856,
                                                                       21812, 21912, 38762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 55082, 0, 3,
                                                                       52856, 37322, 52982,
                                                                       21912, 22012, 38912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 55292, 0, 3,
                                                                       52982, 37412, 53108,
                                                                       22012, 22112, 39062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 55502, 0, 3,
                                                                       53108, 37502, 53234,
                                                                       22112, 22212, 39212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 55712, 0, 3,
                                                                       53234, 37592, 53360,
                                                                       22212, 22312, 39362,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 55922, 0, 3,
                                                                       53360, 37682, 53486,
                                                                       22312, 22412, 39512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 56132, 0, 3,
                                                                       53612, 37862, 53822,
                                                                       22612, 22762, 39662,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 56447, 0, 3,
                                                                       53822, 38012, 54032,
                                                                       22762, 22912, 39887,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 56762, 0, 3,
                                                                       54032, 38162, 54242,
                                                                       22912, 23062, 40112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 57077, 0, 3,
                                                                       54242, 38312, 54452,
                                                                       23062, 23212, 40337,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 57392, 0, 3,
                                                                       54452, 38462, 54662,
                                                                       23212, 23362, 40562,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 57707, 0, 3,
                                                                       54872, 38762, 55082,
                                                                       23662, 23812, 40787,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 58022, 0, 3,
                                                                       55082, 38912, 55292,
                                                                       23812, 23962, 41012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 58337, 0, 3,
                                                                       55292, 39062, 55502,
                                                                       23962, 24112, 41237,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 58652, 0, 3,
                                                                       55502, 39212, 55712,
                                                                       24112, 24262, 41462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 58967, 0, 3,
                                                                       55712, 39362, 55922,
                                                                       24262, 24412, 41687,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 59282, 0, 3,
                                                                       56132, 39662, 56447,
                                                                       24712, 24922, 41912,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 59723, 0, 3,
                                                                       56447, 39887, 56762,
                                                                       24922, 25132, 42227,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 60164, 0, 3,
                                                                       56762, 40112, 57077,
                                                                       25132, 25342, 42542,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 60605, 0, 3,
                                                                       57077, 40337, 57392,
                                                                       25342, 25552, 42857,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 61046, 0, 3,
                                                                       57707, 40787, 58022,
                                                                       25972, 26182, 43172,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 61487, 0, 3,
                                                                       58022, 41012, 58337,
                                                                       26182, 26392, 43487,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 61928, 0, 3,
                                                                       58337, 41237, 58652,
                                                                       26392, 26602, 43802,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 62369, 0, 3,
                                                                       58652, 41462, 58967,
                                                                       26602, 26812, 44117,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 62810, 0, 3,
                                                                       59282, 41912, 59723,
                                                                       27232, 27512, 44432,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 63398, 0, 3,
                                                                       59723, 42227, 60164,
                                                                       27512, 27792, 44852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 63986, 0, 3,
                                                                       60164, 42542, 60605,
                                                                       27792, 28072, 45272,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 64574, 0, 3,
                                                                       61046, 43172, 61487,
                                                                       28632, 28912, 45692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 65162, 0, 3,
                                                                       61487, 43487, 61928,
                                                                       28912, 29192, 46112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 65750, 0, 3,
                                                                       61928, 43802, 62369,
                                                                       29192, 29472, 46532,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 66338, 0, 3,
                                                                       62810, 44432, 63398,
                                                                       30032, 30392, 46952,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 67094, 0, 3,
                                                                       63398, 44852, 63986,
                                                                       30392, 30752, 47492,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 67850, 0, 3,
                                                                       64574, 45692, 65162,
                                                                       31472, 31832, 48032,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 68606, 0, 3,
                                                                       65162, 46112, 65750,
                                                                       31832, 32192, 48572,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 69362, 0, 3,
                                                                       66338, 46952, 67094,
                                                                       32912, 33362, 49112,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 70307, 0, 3,
                                                                       67850, 48032, 68606,
                                                                       34262, 34712, 49787,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_f_x(buffer, 71252, 51848, 56132, 1, 21, ncols, beta);

                    simdgeo::geom_f_y(buffer, 71462, 51848, 56132, 1, 21, ncols, beta);

                    simdgeo::geom_f_z(buffer, 71672, 51848, 56132, 1, 21, ncols, beta);

                    simdgeo::geom_f_x(buffer, 71882, 52730, 57707, 1, 21, ncols, beta);

                    simdgeo::geom_f_y(buffer, 72092, 52730, 57707, 1, 21, ncols, beta);

                    simdgeo::geom_f_z(buffer, 72302, 52730, 57707, 1, 21, ncols, beta);

                    simdgeo::geom_g_x(buffer, 72512, 53612, 59282, 1, 21, ncols, beta);

                    simdgeo::geom_g_y(buffer, 72827, 53612, 59282, 1, 21, ncols, beta);

                    simdgeo::geom_g_z(buffer, 73142, 53612, 59282, 1, 21, ncols, beta);

                    simdgeo::geom_g_x(buffer, 73457, 54872, 61046, 1, 21, ncols, beta);

                    simdgeo::geom_g_y(buffer, 73772, 54872, 61046, 1, 21, ncols, beta);

                    simdgeo::geom_g_z(buffer, 74087, 54872, 61046, 1, 21, ncols, beta);

                    simdgeo::geom_h_x(buffer, 74402, 56132, 62810, 1, 21, ncols, beta);

                    simdgeo::geom_h_y(buffer, 74843, 56132, 62810, 1, 21, ncols, beta);

                    simdgeo::geom_h_z(buffer, 75284, 56132, 62810, 1, 21, ncols, beta);

                    simdgeo::geom_h_x(buffer, 75725, 57707, 64574, 1, 21, ncols, beta);

                    simdgeo::geom_h_y(buffer, 76166, 57707, 64574, 1, 21, ncols, beta);

                    simdgeo::geom_h_z(buffer, 76607, 57707, 64574, 1, 21, ncols, beta);

                    simdgeo::geom_i_x(buffer, 77048, 59282, 66338, 1, 21, ncols, beta);

                    simdgeo::geom_i_y(buffer, 77636, 59282, 66338, 1, 21, ncols, beta);

                    simdgeo::geom_i_z(buffer, 78224, 59282, 66338, 1, 21, ncols, beta);

                    simdgeo::geom_i_x(buffer, 78812, 61046, 67850, 1, 21, ncols, beta);

                    simdgeo::geom_i_y(buffer, 79400, 61046, 67850, 1, 21, ncols, beta);

                    simdgeo::geom_i_z(buffer, 79988, 61046, 67850, 1, 21, ncols, beta);

                    simdgeo::geom_k_x(buffer, 80576, 62810, 69362, 1, 21, ncols, beta);

                    simdgeo::geom_k_y(buffer, 81332, 62810, 69362, 1, 21, ncols, beta);

                    simdgeo::geom_k_z(buffer, 82088, 62810, 69362, 1, 21, ncols, beta);

                    simdgeo::geom_k_x(buffer, 82844, 64574, 70307, 1, 21, ncols, beta);

                    simdgeo::geom_k_y(buffer, 83600, 64574, 70307, 1, 21, ncols, beta);

                    simdgeo::geom_k_z(buffer, 84356, 64574, 70307, 1, 21, ncols, beta);

                    simdfunc::contract_primitives(buffer, 85112, 71252, 210, ncols);

                    simdfunc::contract_primitives(buffer, 85432, 71462, 210, ncols);

                    simdfunc::contract_primitives(buffer, 85752, 71672, 210, ncols);

                    simdfunc::contract_primitives(buffer, 86072, 53612, 210, ncols);

                    simdfunc::contract_primitives(buffer, 86392, 71882, 210, ncols);

                    simdfunc::contract_primitives(buffer, 86712, 72092, 210, ncols);

                    simdfunc::contract_primitives(buffer, 87032, 72302, 210, ncols);

                    simdfunc::contract_primitives(buffer, 87352, 54872, 210, ncols);

                    simdfunc::contract_primitives(buffer, 87672, 72512, 315, ncols);

                    simdfunc::contract_primitives(buffer, 88152, 72827, 315, ncols);

                    simdfunc::contract_primitives(buffer, 88632, 73142, 315, ncols);

                    simdfunc::contract_primitives(buffer, 89112, 56132, 315, ncols);

                    simdfunc::contract_primitives(buffer, 89592, 73457, 315, ncols);

                    simdfunc::contract_primitives(buffer, 90072, 73772, 315, ncols);

                    simdfunc::contract_primitives(buffer, 90552, 74087, 315, ncols);

                    simdfunc::contract_primitives(buffer, 91032, 57707, 315, ncols);

                    simdfunc::contract_primitives(buffer, 91512, 74402, 441, ncols);

                    simdfunc::contract_primitives(buffer, 92184, 74843, 441, ncols);

                    simdfunc::contract_primitives(buffer, 92856, 75284, 441, ncols);

                    simdfunc::contract_primitives(buffer, 93528, 59282, 441, ncols);

                    simdfunc::contract_primitives(buffer, 94200, 75725, 441, ncols);

                    simdfunc::contract_primitives(buffer, 94872, 76166, 441, ncols);

                    simdfunc::contract_primitives(buffer, 95544, 76607, 441, ncols);

                    simdfunc::contract_primitives(buffer, 96216, 61046, 441, ncols);

                    simdfunc::contract_primitives(buffer, 96888, 77048, 588, ncols);

                    simdfunc::contract_primitives(buffer, 97784, 77636, 588, ncols);

                    simdfunc::contract_primitives(buffer, 98680, 78224, 588, ncols);

                    simdfunc::contract_primitives(buffer, 99576, 62810, 588, ncols);

                    simdfunc::contract_primitives(buffer, 100472, 78812, 588, ncols);

                    simdfunc::contract_primitives(buffer, 101368, 79400, 588, ncols);

                    simdfunc::contract_primitives(buffer, 102264, 79988, 588, ncols);

                    simdfunc::contract_primitives(buffer, 103160, 64574, 588, ncols);

                    simdfunc::contract_primitives(buffer, 104056, 80576, 756, ncols);

                    simdfunc::contract_primitives(buffer, 105208, 81332, 756, ncols);

                    simdfunc::contract_primitives(buffer, 106360, 82088, 756, ncols);

                    simdfunc::contract_primitives(buffer, 107512, 82844, 756, ncols);

                    simdfunc::contract_primitives(buffer, 108664, 83600, 756, ncols);

                    simdfunc::contract_primitives(buffer, 109816, 84356, 756, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 85322, 85112, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 85642, 85432, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 85962, 85752, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 86282, 86072, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 86602, 86392, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 86922, 86712, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 87242, 87032, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 87562, 87352, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 87987, 87672, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 88467, 88152, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 88947, 88632, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 89427, 89112, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 89907, 89592, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 90387, 90072, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 90867, 90552, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 91347, 91032, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 91953, 91512, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 92625, 92184, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 93297, 92856, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 93969, 93528, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 94641, 94200, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 95313, 94872, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 95985, 95544, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 96657, 96216, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 97476, 96888, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 98372, 97784, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 99268, 98680, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 100164, 99576, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 101060, 100472, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 101956, 101368, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 102852, 102264, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 103748, 103160, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 104812, 104056, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 105964, 105208, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 107116, 106360, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 108268, 107512, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 109420, 108664, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 110572, 109816, 36, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 110968, 85322, 86282, 87987, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 111298, 85642, 86282, 88467, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 111628, 85962, 86282, 88947, 11,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 111958, 86282, 89427, 11, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 112288, 86602, 87562, 89907, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 112618, 86922, 87562, 90387, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 112948, 87242, 87562, 90867, 11,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 113278, 87562, 91347, 11, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 113608, 87987, 89427, 91953, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 114103, 88467, 89427, 92625, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 114598, 88947, 89427, 93297, 11,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 115093, 89427, 93969, 11, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 115588, 89907, 91347, 94641, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 116083, 90387, 91347, 95313, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 116578, 90867, 91347, 95985, 11,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 117073, 91347, 96657, 11, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 117568, 91953, 93969, 97476, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 118261, 92625, 93969, 98372, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 118954, 93297, 93969, 99268, 11,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 119647, 93969, 100164, 11, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 120340, 94641, 96657, 101060, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 121033, 95313, 96657, 101956, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 121726, 95985, 96657, 102852, 11,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 122419, 96657, 103748, 11, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 123112, 97476, 100164, 104812, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 124036, 98372, 100164, 105964, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 124960, 99268, 100164, 107116, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 125884, 101060, 103748, 108268,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 126808, 101956, 103748, 109420,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 127732, 102852, 103748, 110572,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 128656, 110968, 111958, 113608,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 129316, 111298, 111958, 114103,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 129976, 111628, 111958, 114598,
                                          11, nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 130636, 111958, 115093, 11, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 131296, 112288, 113278, 115588,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 131956, 112618, 113278, 116083,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 132616, 112948, 113278, 116578,
                                          11, nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 133276, 113278, 117073, 11, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 133936, 113608, 115093, 117568,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 134926, 114103, 115093, 118261,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 135916, 114598, 115093, 118954,
                                          11, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 136906, 115093, 119647, 11, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 137896, 115588, 117073, 120340,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 138886, 116083, 117073, 121033,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 139876, 116578, 117073, 121726,
                                          11, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 140866, 117073, 122419, 11, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 141856, 117568, 119647, 123112,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 143242, 118261, 119647, 124036,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 144628, 118954, 119647, 124960,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 146014, 120340, 122419, 125884,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 147400, 121033, 122419, 126808,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 148786, 121726, 122419, 127732,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 150172, 128656, 130636, 133936,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 151272, 129316, 130636, 134926,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 152372, 129976, 130636, 135916,
                                          11, nmax);

        simdtrf::compute_hrr_ff(buffer, coordinates, 153472, 130636, 136906, 11, nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 154572, 131296, 133276, 137896,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 155672, 131956, 133276, 138886,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 156772, 132616, 133276, 139876,
                                          11, nmax);

        simdtrf::compute_hrr_ff(buffer, coordinates, 157872, 133276, 140866, 11, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 158972, 133936, 136906, 141856,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 160622, 134926, 136906, 143242,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 162272, 135916, 136906, 144628,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 163922, 137896, 140866, 146014,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 165572, 138886, 140866, 147400,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 167222, 139876, 140866, 148786,
                                          11, nmax);

        simdtrf::compute_hrr_geom_010x_gf_out_of_second(buffer, coordinates, 168872, 150172,
                                                        153472, 158972, 11, nmax);

        simdtrf::compute_hrr_geom_010y_gf_out_of_second(buffer, coordinates, 170522, 151272,
                                                        153472, 160622, 11, nmax);

        simdtrf::compute_hrr_geom_010z_gf_out_of_second(buffer, coordinates, 172172, 152372,
                                                        153472, 162272, 11, nmax);

        simdtrf::compute_hrr_geom_010x_gf_out_of_second(buffer, coordinates, 173822, 154572,
                                                        157872, 163922, 11, nmax);

        simdtrf::compute_hrr_geom_010y_gf_out_of_second(buffer, coordinates, 175472, 155672,
                                                        157872, 165572, 11, nmax);

        simdtrf::compute_hrr_geom_010z_gf_out_of_second(buffer, coordinates, 177122, 156772,
                                                        157872, 167222, 11, nmax);

        simdtrf::transform_f_inner(buffer, 178772, 173822, 15, 11, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 178772, 77, nmax);

        simdtrf::transform_f_inner(buffer, 178772, 175472, 15, 11, nmax);

        simdtrf::transform_g_outer(values + 693 * nvalues + n * npairs, nvalues, buffer, 178772,
                                   77, nmax);

        simdtrf::transform_f_inner(buffer, 178772, 177122, 15, 11, nmax);

        simdtrf::transform_g_outer(values + 1386 * nvalues + n * npairs, nvalues, buffer, 178772,
                                   77, nmax);

        simdtrf::transform_f_inner(buffer, 178772, 168872, 15, 11, nmax);

        simdtrf::transform_g_outer(values + 2079 * nvalues + n * npairs, nvalues, buffer, 178772,
                                   77, nmax);

        simdtrf::transform_f_inner(buffer, 178772, 170522, 15, 11, nmax);

        simdtrf::transform_g_outer(values + 2772 * nvalues + n * npairs, nvalues, buffer, 178772,
                                   77, nmax);

        simdtrf::transform_f_inner(buffer, 178772, 172172, 15, 11, nmax);

        simdtrf::transform_g_outer(values + 3465 * nvalues + n * npairs, nvalues, buffer, 178772,
                                   77, nmax);
    }

    for (size_t m = 0; m < 4158; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
