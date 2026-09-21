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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGFG.hpp"

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
compute_rs_geom_010_gfg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_gfg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 127631, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 3402 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 127631, 51818, 19068, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 12,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 20, 3, 12,
                                                             ncols, fj, i * nprim_b + j, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2318, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2321, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2324, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2327, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2330, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2333, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2336, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2339, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2342, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2345, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2348, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2351, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2354, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2357, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2360, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2363, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2366, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2369, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2372, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2375, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2378, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2381, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2384, 3, 9, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2393, 3, 10, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2402, 3, 11, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2411, 3, 12, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2420, 3, 13, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2429, 3, 14, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2438, 3, 15, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2447, 3, 16, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2456, 3, 17, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2465, 3, 18, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2474, 3, 23, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2483, 3, 24, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2492, 3, 25, 82,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2501, 3, 26, 85,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2510, 3, 27, 88,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2519, 3, 28, 91,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2528, 3, 29, 94,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2537, 3, 30, 97,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2546, 3, 31, 100,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2555, 3, 32, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2564, 3, 40, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2582, 3, 43, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2600, 3, 46, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2618, 3, 49, 136,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2636, 3, 52, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2654, 3, 55, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2672, 3, 58, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2690, 3, 61, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2708, 3, 64, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2726, 3, 76, 184,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2744, 3, 79, 190,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2762, 3, 82, 196,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2780, 3, 85, 202,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2798, 3, 88, 208,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2816, 3, 91, 214,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2834, 3, 94, 220,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2852, 3, 97, 226,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2870, 3, 100, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2888, 3, 118, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2918, 3, 124, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2948, 3, 130, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2978, 3, 136, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3008, 3, 142, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3038, 3, 148, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3068, 3, 154, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3098, 3, 160, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3128, 3, 184, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3158, 3, 190, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3188, 3, 196, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3218, 3, 202, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3248, 3, 208, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3278, 3, 214, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3308, 3, 220, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3338, 3, 226, 428,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3368, 3, 258, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3413, 3, 268, 483,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3458, 3, 278, 498,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3503, 3, 288, 513,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3548, 3, 298, 528,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3593, 3, 308, 543,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3638, 3, 318, 558,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3683, 3, 358, 603,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3728, 3, 368, 618,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3773, 3, 378, 633,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3818, 3, 388, 648,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3863, 3, 398, 663,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3908, 3, 408, 678,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3953, 3, 418, 693,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3998, 3, 468, 750,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4061, 3, 483, 771,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4124, 3, 498, 792,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4187, 3, 513, 813,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4250, 3, 528, 834,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4313, 3, 543, 855,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4376, 3, 603, 918,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4439, 3, 618, 939,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4502, 3, 633, 960,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4565, 3, 648, 981,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4628, 3, 663,
                                                                       1002, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4691, 3, 678,
                                                                       1023, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4754, 3, 750,
                                                                       1100, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4838, 3, 771,
                                                                       1128, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4922, 3, 792,
                                                                       1156, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5006, 3, 813,
                                                                       1184, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5090, 3, 834,
                                                                       1212, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5174, 3, 918,
                                                                       1296, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5258, 3, 939,
                                                                       1324, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5342, 3, 960,
                                                                       1352, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5426, 3, 981,
                                                                       1380, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5510, 3, 1002,
                                                                       1408, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5594, 3, 1100,
                                                                       1508, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5702, 3, 1128,
                                                                       1544, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5810, 3, 1156,
                                                                       1580, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5918, 3, 1184,
                                                                       1616, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6026, 3, 1296,
                                                                       1724, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6134, 3, 1324,
                                                                       1760, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6242, 3, 1352,
                                                                       1796, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6350, 3, 1380,
                                                                       1832, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6458, 3, 1508,
                                                                       1958, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6593, 3, 1544,
                                                                       2003, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6728, 3, 1580,
                                                                       2048, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6863, 3, 1724,
                                                                       2183, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6998, 3, 1760,
                                                                       2228, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7133, 3, 1796,
                                                                       2273, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7268, 3, 7, 8,
                                                                       2318, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7274, 3, 8, 9,
                                                                       2321, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7280, 3, 9, 10,
                                                                       2324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7286, 3, 10, 11,
                                                                       2327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7292, 3, 11, 12,
                                                                       2330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7298, 3, 12, 13,
                                                                       2333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7304, 3, 13, 14,
                                                                       2336, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7310, 3, 14, 15,
                                                                       2339, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7316, 3, 15, 16,
                                                                       2342, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7322, 3, 16, 17,
                                                                       2345, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7328, 3, 17, 18,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7334, 3, 21, 22,
                                                                       2351, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7340, 3, 22, 23,
                                                                       2354, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7346, 3, 23, 24,
                                                                       2357, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7352, 3, 24, 25,
                                                                       2360, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7358, 3, 25, 26,
                                                                       2363, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7364, 3, 26, 27,
                                                                       2366, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7370, 3, 27, 28,
                                                                       2369, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7376, 3, 28, 29,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7382, 3, 29, 30,
                                                                       2375, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7388, 3, 30, 31,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7394, 3, 31, 32,
                                                                       2381, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 7268,
                                                                       2318, 7274, 2384, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7418, 0, 3, 7274,
                                                                       2321, 7280, 2393, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7436, 0, 3, 7280,
                                                                       2324, 7286, 2402, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7454, 0, 3, 7286,
                                                                       2327, 7292, 2411, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7472, 0, 3, 7292,
                                                                       2330, 7298, 2420, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7490, 0, 3, 7298,
                                                                       2333, 7304, 2429, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7508, 0, 3, 7304,
                                                                       2336, 7310, 2438, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7526, 0, 3, 7310,
                                                                       2339, 7316, 2447, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7544, 0, 3, 7316,
                                                                       2342, 7322, 2456, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7562, 0, 3, 7322,
                                                                       2345, 7328, 2465, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7580, 0, 3, 7334,
                                                                       2351, 7340, 2474, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7598, 0, 3, 7340,
                                                                       2354, 7346, 2483, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7616, 0, 3, 7346,
                                                                       2357, 7352, 2492, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7634, 0, 3, 7352,
                                                                       2360, 7358, 2501, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7652, 0, 3, 7358,
                                                                       2363, 7364, 2510, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7670, 0, 3, 7364,
                                                                       2366, 7370, 2519, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7688, 0, 3, 7370,
                                                                       2369, 7376, 2528, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7706, 0, 3, 7376,
                                                                       2372, 7382, 2537, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7724, 0, 3, 7382,
                                                                       2375, 7388, 2546, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7742, 0, 3, 7388,
                                                                       2378, 7394, 2555, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7760, 0, 3, 7400,
                                                                       2384, 7418, 106, 112,
                                                                       2564, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7796, 0, 3, 7418,
                                                                       2393, 7436, 112, 118,
                                                                       2582, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7832, 0, 3, 7436,
                                                                       2402, 7454, 118, 124,
                                                                       2600, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7868, 0, 3, 7454,
                                                                       2411, 7472, 124, 130,
                                                                       2618, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7904, 0, 3, 7472,
                                                                       2420, 7490, 130, 136,
                                                                       2636, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7940, 0, 3, 7490,
                                                                       2429, 7508, 136, 142,
                                                                       2654, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7976, 0, 3, 7508,
                                                                       2438, 7526, 142, 148,
                                                                       2672, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8012, 0, 3, 7526,
                                                                       2447, 7544, 148, 154,
                                                                       2690, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8048, 0, 3, 7544,
                                                                       2456, 7562, 154, 160,
                                                                       2708, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8084, 0, 3, 7580,
                                                                       2474, 7598, 172, 178,
                                                                       2726, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8120, 0, 3, 7598,
                                                                       2483, 7616, 178, 184,
                                                                       2744, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8156, 0, 3, 7616,
                                                                       2492, 7634, 184, 190,
                                                                       2762, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8192, 0, 3, 7634,
                                                                       2501, 7652, 190, 196,
                                                                       2780, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8228, 0, 3, 7652,
                                                                       2510, 7670, 196, 202,
                                                                       2798, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8264, 0, 3, 7670,
                                                                       2519, 7688, 202, 208,
                                                                       2816, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8300, 0, 3, 7688,
                                                                       2528, 7706, 208, 214,
                                                                       2834, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8336, 0, 3, 7706,
                                                                       2537, 7724, 214, 220,
                                                                       2852, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8372, 0, 3, 7724,
                                                                       2546, 7742, 220, 226,
                                                                       2870, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8408, 0, 3, 7760,
                                                                       2564, 7796, 238, 248,
                                                                       2888, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8468, 0, 3, 7796,
                                                                       2582, 7832, 248, 258,
                                                                       2918, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8528, 0, 3, 7832,
                                                                       2600, 7868, 258, 268,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8588, 0, 3, 7868,
                                                                       2618, 7904, 268, 278,
                                                                       2978, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8648, 0, 3, 7904,
                                                                       2636, 7940, 278, 288,
                                                                       3008, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8708, 0, 3, 7940,
                                                                       2654, 7976, 288, 298,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8768, 0, 3, 7976,
                                                                       2672, 8012, 298, 308,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8828, 0, 3, 8012,
                                                                       2690, 8048, 308, 318,
                                                                       3098, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8888, 0, 3, 8084,
                                                                       2726, 8120, 338, 348,
                                                                       3128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8948, 0, 3, 8120,
                                                                       2744, 8156, 348, 358,
                                                                       3158, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9008, 0, 3, 8156,
                                                                       2762, 8192, 358, 368,
                                                                       3188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9068, 0, 3, 8192,
                                                                       2780, 8228, 368, 378,
                                                                       3218, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9128, 0, 3, 8228,
                                                                       2798, 8264, 378, 388,
                                                                       3248, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9188, 0, 3, 8264,
                                                                       2816, 8300, 388, 398,
                                                                       3278, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9248, 0, 3, 8300,
                                                                       2834, 8336, 398, 408,
                                                                       3308, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9308, 0, 3, 8336,
                                                                       2852, 8372, 408, 418,
                                                                       3338, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9368, 0, 3, 8408,
                                                                       2888, 8468, 438, 453,
                                                                       3368, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9458, 0, 3, 8468,
                                                                       2918, 8528, 453, 468,
                                                                       3413, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9548, 0, 3, 8528,
                                                                       2948, 8588, 468, 483,
                                                                       3458, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9638, 0, 3, 8588,
                                                                       2978, 8648, 483, 498,
                                                                       3503, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9728, 0, 3, 8648,
                                                                       3008, 8708, 498, 513,
                                                                       3548, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9818, 0, 3, 8708,
                                                                       3038, 8768, 513, 528,
                                                                       3593, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9908, 0, 3, 8768,
                                                                       3068, 8828, 528, 543,
                                                                       3638, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9998, 0, 3, 8888,
                                                                       3128, 8948, 573, 588,
                                                                       3683, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10088, 0, 3, 8948,
                                                                       3158, 9008, 588, 603,
                                                                       3728, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10178, 0, 3, 9008,
                                                                       3188, 9068, 603, 618,
                                                                       3773, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10268, 0, 3, 9068,
                                                                       3218, 9128, 618, 633,
                                                                       3818, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10358, 0, 3, 9128,
                                                                       3248, 9188, 633, 648,
                                                                       3863, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10448, 0, 3, 9188,
                                                                       3278, 9248, 648, 663,
                                                                       3908, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10538, 0, 3, 9248,
                                                                       3308, 9308, 663, 678,
                                                                       3953, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10628, 0, 3, 9368,
                                                                       3368, 9458, 708, 729,
                                                                       3998, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10754, 0, 3, 9458,
                                                                       3413, 9548, 729, 750,
                                                                       4061, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10880, 0, 3, 9548,
                                                                       3458, 9638, 750, 771,
                                                                       4124, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11006, 0, 3, 9638,
                                                                       3503, 9728, 771, 792,
                                                                       4187, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11132, 0, 3, 9728,
                                                                       3548, 9818, 792, 813,
                                                                       4250, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11258, 0, 3, 9818,
                                                                       3593, 9908, 813, 834,
                                                                       4313, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11384, 0, 3, 9998,
                                                                       3683, 10088, 876, 897,
                                                                       4376, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11510, 0, 3,
                                                                       10088, 3728, 10178, 897,
                                                                       918, 4439, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11636, 0, 3,
                                                                       10178, 3773, 10268, 918,
                                                                       939, 4502, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11762, 0, 3,
                                                                       10268, 3818, 10358, 939,
                                                                       960, 4565, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11888, 0, 3,
                                                                       10358, 3863, 10448, 960,
                                                                       981, 4628, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12014, 0, 3,
                                                                       10448, 3908, 10538, 981,
                                                                       1002, 4691, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12140, 0, 3,
                                                                       10628, 3998, 10754, 1044,
                                                                       1072, 4754, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12308, 0, 3,
                                                                       10754, 4061, 10880, 1072,
                                                                       1100, 4838, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12476, 0, 3,
                                                                       10880, 4124, 11006, 1100,
                                                                       1128, 4922, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12644, 0, 3,
                                                                       11006, 4187, 11132, 1128,
                                                                       1156, 5006, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12812, 0, 3,
                                                                       11132, 4250, 11258, 1156,
                                                                       1184, 5090, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12980, 0, 3,
                                                                       11384, 4376, 11510, 1240,
                                                                       1268, 5174, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13148, 0, 3,
                                                                       11510, 4439, 11636, 1268,
                                                                       1296, 5258, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13316, 0, 3,
                                                                       11636, 4502, 11762, 1296,
                                                                       1324, 5342, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13484, 0, 3,
                                                                       11762, 4565, 11888, 1324,
                                                                       1352, 5426, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13652, 0, 3,
                                                                       11888, 4628, 12014, 1352,
                                                                       1380, 5510, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13820, 0, 3,
                                                                       12140, 4754, 12308, 1436,
                                                                       1472, 5594, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14036, 0, 3,
                                                                       12308, 4838, 12476, 1472,
                                                                       1508, 5702, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14252, 0, 3,
                                                                       12476, 4922, 12644, 1508,
                                                                       1544, 5810, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14468, 0, 3,
                                                                       12644, 5006, 12812, 1544,
                                                                       1580, 5918, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14684, 0, 3,
                                                                       12980, 5174, 13148, 1652,
                                                                       1688, 6026, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14900, 0, 3,
                                                                       13148, 5258, 13316, 1688,
                                                                       1724, 6134, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15116, 0, 3,
                                                                       13316, 5342, 13484, 1724,
                                                                       1760, 6242, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15332, 0, 3,
                                                                       13484, 5426, 13652, 1760,
                                                                       1796, 6350, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 15548, 0, 3,
                                                                       13820, 5594, 14036, 1868,
                                                                       1913, 6458, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 15818, 0, 3,
                                                                       14036, 5702, 14252, 1913,
                                                                       1958, 6593, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16088, 0, 3,
                                                                       14252, 5810, 14468, 1958,
                                                                       2003, 6728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16358, 0, 3,
                                                                       14684, 6026, 14900, 2093,
                                                                       2138, 6863, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16628, 0, 3,
                                                                       14900, 6134, 15116, 2138,
                                                                       2183, 6998, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16898, 0, 3,
                                                                       15116, 6242, 15332, 2183,
                                                                       2228, 7133, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17168, 3, 2318,
                                                                       2321, 7280, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17178, 3, 2321,
                                                                       2324, 7286, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17188, 3, 2324,
                                                                       2327, 7292, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17198, 3, 2327,
                                                                       2330, 7298, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17208, 3, 2330,
                                                                       2333, 7304, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17218, 3, 2333,
                                                                       2336, 7310, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17228, 3, 2336,
                                                                       2339, 7316, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17238, 3, 2339,
                                                                       2342, 7322, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17248, 3, 2342,
                                                                       2345, 7328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17258, 3, 2351,
                                                                       2354, 7346, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17268, 3, 2354,
                                                                       2357, 7352, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17278, 3, 2357,
                                                                       2360, 7358, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17288, 3, 2360,
                                                                       2363, 7364, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17298, 3, 2363,
                                                                       2366, 7370, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17308, 3, 2366,
                                                                       2369, 7376, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17318, 3, 2369,
                                                                       2372, 7382, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17328, 3, 2372,
                                                                       2375, 7388, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17338, 3, 2375,
                                                                       2378, 7394, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17348, 0, 3,
                                                                       17168, 7280, 17178, 7436,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17378, 0, 3,
                                                                       17178, 7286, 17188, 7454,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17408, 0, 3,
                                                                       17188, 7292, 17198, 7472,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17438, 0, 3,
                                                                       17198, 7298, 17208, 7490,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17468, 0, 3,
                                                                       17208, 7304, 17218, 7508,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17498, 0, 3,
                                                                       17218, 7310, 17228, 7526,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17528, 0, 3,
                                                                       17228, 7316, 17238, 7544,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17558, 0, 3,
                                                                       17238, 7322, 17248, 7562,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17588, 0, 3,
                                                                       17258, 7346, 17268, 7616,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17618, 0, 3,
                                                                       17268, 7352, 17278, 7634,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17648, 0, 3,
                                                                       17278, 7358, 17288, 7652,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17678, 0, 3,
                                                                       17288, 7364, 17298, 7670,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17708, 0, 3,
                                                                       17298, 7370, 17308, 7688,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17738, 0, 3,
                                                                       17308, 7376, 17318, 7706,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17768, 0, 3,
                                                                       17318, 7382, 17328, 7724,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17798, 0, 3,
                                                                       17328, 7388, 17338, 7742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17828, 0, 3,
                                                                       17348, 7436, 17378, 2564,
                                                                       2582, 7832, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17888, 0, 3,
                                                                       17378, 7454, 17408, 2582,
                                                                       2600, 7868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17948, 0, 3,
                                                                       17408, 7472, 17438, 2600,
                                                                       2618, 7904, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18008, 0, 3,
                                                                       17438, 7490, 17468, 2618,
                                                                       2636, 7940, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18068, 0, 3,
                                                                       17468, 7508, 17498, 2636,
                                                                       2654, 7976, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18128, 0, 3,
                                                                       17498, 7526, 17528, 2654,
                                                                       2672, 8012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18188, 0, 3,
                                                                       17528, 7544, 17558, 2672,
                                                                       2690, 8048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18248, 0, 3,
                                                                       17588, 7616, 17618, 2726,
                                                                       2744, 8156, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18308, 0, 3,
                                                                       17618, 7634, 17648, 2744,
                                                                       2762, 8192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18368, 0, 3,
                                                                       17648, 7652, 17678, 2762,
                                                                       2780, 8228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18428, 0, 3,
                                                                       17678, 7670, 17708, 2780,
                                                                       2798, 8264, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18488, 0, 3,
                                                                       17708, 7688, 17738, 2798,
                                                                       2816, 8300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18548, 0, 3,
                                                                       17738, 7706, 17768, 2816,
                                                                       2834, 8336, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18608, 0, 3,
                                                                       17768, 7724, 17798, 2834,
                                                                       2852, 8372, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18668, 0, 3,
                                                                       17828, 7832, 17888, 2888,
                                                                       2918, 8528, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18768, 0, 3,
                                                                       17888, 7868, 17948, 2918,
                                                                       2948, 8588, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18868, 0, 3,
                                                                       17948, 7904, 18008, 2948,
                                                                       2978, 8648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18968, 0, 3,
                                                                       18008, 7940, 18068, 2978,
                                                                       3008, 8708, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19068, 0, 3,
                                                                       18068, 7976, 18128, 3008,
                                                                       3038, 8768, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19168, 0, 3,
                                                                       18128, 8012, 18188, 3038,
                                                                       3068, 8828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19268, 0, 3,
                                                                       18248, 8156, 18308, 3128,
                                                                       3158, 9008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19368, 0, 3,
                                                                       18308, 8192, 18368, 3158,
                                                                       3188, 9068, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19468, 0, 3,
                                                                       18368, 8228, 18428, 3188,
                                                                       3218, 9128, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19568, 0, 3,
                                                                       18428, 8264, 18488, 3218,
                                                                       3248, 9188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19668, 0, 3,
                                                                       18488, 8300, 18548, 3248,
                                                                       3278, 9248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19768, 0, 3,
                                                                       18548, 8336, 18608, 3278,
                                                                       3308, 9308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19868, 0, 3,
                                                                       18668, 8528, 18768, 3368,
                                                                       3413, 9548, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20018, 0, 3,
                                                                       18768, 8588, 18868, 3413,
                                                                       3458, 9638, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20168, 0, 3,
                                                                       18868, 8648, 18968, 3458,
                                                                       3503, 9728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20318, 0, 3,
                                                                       18968, 8708, 19068, 3503,
                                                                       3548, 9818, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20468, 0, 3,
                                                                       19068, 8768, 19168, 3548,
                                                                       3593, 9908, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20618, 0, 3,
                                                                       19268, 9008, 19368, 3683,
                                                                       3728, 10178, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20768, 0, 3,
                                                                       19368, 9068, 19468, 3728,
                                                                       3773, 10268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20918, 0, 3,
                                                                       19468, 9128, 19568, 3773,
                                                                       3818, 10358, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21068, 0, 3,
                                                                       19568, 9188, 19668, 3818,
                                                                       3863, 10448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21218, 0, 3,
                                                                       19668, 9248, 19768, 3863,
                                                                       3908, 10538, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21368, 0, 3,
                                                                       19868, 9548, 20018, 3998,
                                                                       4061, 10880, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21578, 0, 3,
                                                                       20018, 9638, 20168, 4061,
                                                                       4124, 11006, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21788, 0, 3,
                                                                       20168, 9728, 20318, 4124,
                                                                       4187, 11132, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21998, 0, 3,
                                                                       20318, 9818, 20468, 4187,
                                                                       4250, 11258, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22208, 0, 3,
                                                                       20618, 10178, 20768, 4376,
                                                                       4439, 11636, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22418, 0, 3,
                                                                       20768, 10268, 20918, 4439,
                                                                       4502, 11762, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22628, 0, 3,
                                                                       20918, 10358, 21068, 4502,
                                                                       4565, 11888, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22838, 0, 3,
                                                                       21068, 10448, 21218, 4565,
                                                                       4628, 12014, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23048, 0, 3,
                                                                       21368, 10880, 21578, 4754,
                                                                       4838, 12476, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23328, 0, 3,
                                                                       21578, 11006, 21788, 4838,
                                                                       4922, 12644, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23608, 0, 3,
                                                                       21788, 11132, 21998, 4922,
                                                                       5006, 12812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23888, 0, 3,
                                                                       22208, 11636, 22418, 5174,
                                                                       5258, 13316, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24168, 0, 3,
                                                                       22418, 11762, 22628, 5258,
                                                                       5342, 13484, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24448, 0, 3,
                                                                       22628, 11888, 22838, 5342,
                                                                       5426, 13652, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 24728, 0, 3,
                                                                       23048, 12476, 23328, 5594,
                                                                       5702, 14252, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 25088, 0, 3,
                                                                       23328, 12644, 23608, 5702,
                                                                       5810, 14468, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 25448, 0, 3,
                                                                       23888, 13316, 24168, 6026,
                                                                       6134, 15116, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 25808, 0, 3,
                                                                       24168, 13484, 24448, 6134,
                                                                       6242, 15332, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 26168, 0, 3,
                                                                       24728, 14252, 25088, 6458,
                                                                       6593, 16088, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 26618, 0, 3,
                                                                       25448, 15116, 25808, 6863,
                                                                       6998, 16898, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27068, 3, 7268,
                                                                       7274, 17168, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27083, 3, 7274,
                                                                       7280, 17178, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27098, 3, 7280,
                                                                       7286, 17188, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27113, 3, 7286,
                                                                       7292, 17198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27128, 3, 7292,
                                                                       7298, 17208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27143, 3, 7298,
                                                                       7304, 17218, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27158, 3, 7304,
                                                                       7310, 17228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27173, 3, 7310,
                                                                       7316, 17238, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27188, 3, 7316,
                                                                       7322, 17248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27203, 3, 7334,
                                                                       7340, 17258, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27218, 3, 7340,
                                                                       7346, 17268, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27233, 3, 7346,
                                                                       7352, 17278, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27248, 3, 7352,
                                                                       7358, 17288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27263, 3, 7358,
                                                                       7364, 17298, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27278, 3, 7364,
                                                                       7370, 17308, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27293, 3, 7370,
                                                                       7376, 17318, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27308, 3, 7376,
                                                                       7382, 17328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27323, 3, 7382,
                                                                       7388, 17338, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27338, 0, 3,
                                                                       27068, 17168, 27083, 7400,
                                                                       7418, 17348, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27383, 0, 3,
                                                                       27083, 17178, 27098, 7418,
                                                                       7436, 17378, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27428, 0, 3,
                                                                       27098, 17188, 27113, 7436,
                                                                       7454, 17408, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27473, 0, 3,
                                                                       27113, 17198, 27128, 7454,
                                                                       7472, 17438, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27518, 0, 3,
                                                                       27128, 17208, 27143, 7472,
                                                                       7490, 17468, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27563, 0, 3,
                                                                       27143, 17218, 27158, 7490,
                                                                       7508, 17498, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27608, 0, 3,
                                                                       27158, 17228, 27173, 7508,
                                                                       7526, 17528, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27653, 0, 3,
                                                                       27173, 17238, 27188, 7526,
                                                                       7544, 17558, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27698, 0, 3,
                                                                       27203, 17258, 27218, 7580,
                                                                       7598, 17588, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27743, 0, 3,
                                                                       27218, 17268, 27233, 7598,
                                                                       7616, 17618, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27788, 0, 3,
                                                                       27233, 17278, 27248, 7616,
                                                                       7634, 17648, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27833, 0, 3,
                                                                       27248, 17288, 27263, 7634,
                                                                       7652, 17678, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27878, 0, 3,
                                                                       27263, 17298, 27278, 7652,
                                                                       7670, 17708, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27923, 0, 3,
                                                                       27278, 17308, 27293, 7670,
                                                                       7688, 17738, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 27968, 0, 3,
                                                                       27293, 17318, 27308, 7688,
                                                                       7706, 17768, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28013, 0, 3,
                                                                       27308, 17328, 27323, 7706,
                                                                       7724, 17798, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28058, 0, 3,
                                                                       27338, 17348, 27383, 7760,
                                                                       7796, 17828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28148, 0, 3,
                                                                       27383, 17378, 27428, 7796,
                                                                       7832, 17888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28238, 0, 3,
                                                                       27428, 17408, 27473, 7832,
                                                                       7868, 17948, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28328, 0, 3,
                                                                       27473, 17438, 27518, 7868,
                                                                       7904, 18008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28418, 0, 3,
                                                                       27518, 17468, 27563, 7904,
                                                                       7940, 18068, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28508, 0, 3,
                                                                       27563, 17498, 27608, 7940,
                                                                       7976, 18128, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28598, 0, 3,
                                                                       27608, 17528, 27653, 7976,
                                                                       8012, 18188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28688, 0, 3,
                                                                       27698, 17588, 27743, 8084,
                                                                       8120, 18248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28778, 0, 3,
                                                                       27743, 17618, 27788, 8120,
                                                                       8156, 18308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28868, 0, 3,
                                                                       27788, 17648, 27833, 8156,
                                                                       8192, 18368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28958, 0, 3,
                                                                       27833, 17678, 27878, 8192,
                                                                       8228, 18428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 29048, 0, 3,
                                                                       27878, 17708, 27923, 8228,
                                                                       8264, 18488, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 29138, 0, 3,
                                                                       27923, 17738, 27968, 8264,
                                                                       8300, 18548, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 29228, 0, 3,
                                                                       27968, 17768, 28013, 8300,
                                                                       8336, 18608, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29318, 0, 3,
                                                                       28058, 17828, 28148, 8408,
                                                                       8468, 18668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29468, 0, 3,
                                                                       28148, 17888, 28238, 8468,
                                                                       8528, 18768, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29618, 0, 3,
                                                                       28238, 17948, 28328, 8528,
                                                                       8588, 18868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29768, 0, 3,
                                                                       28328, 18008, 28418, 8588,
                                                                       8648, 18968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29918, 0, 3,
                                                                       28418, 18068, 28508, 8648,
                                                                       8708, 19068, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30068, 0, 3,
                                                                       28508, 18128, 28598, 8708,
                                                                       8768, 19168, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30218, 0, 3,
                                                                       28688, 18248, 28778, 8888,
                                                                       8948, 19268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30368, 0, 3,
                                                                       28778, 18308, 28868, 8948,
                                                                       9008, 19368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30518, 0, 3,
                                                                       28868, 18368, 28958, 9008,
                                                                       9068, 19468, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30668, 0, 3,
                                                                       28958, 18428, 29048, 9068,
                                                                       9128, 19568, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30818, 0, 3,
                                                                       29048, 18488, 29138, 9128,
                                                                       9188, 19668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30968, 0, 3,
                                                                       29138, 18548, 29228, 9188,
                                                                       9248, 19768, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31118, 0, 3,
                                                                       29318, 18668, 29468, 9368,
                                                                       9458, 19868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31343, 0, 3,
                                                                       29468, 18768, 29618, 9458,
                                                                       9548, 20018, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31568, 0, 3,
                                                                       29618, 18868, 29768, 9548,
                                                                       9638, 20168, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31793, 0, 3,
                                                                       29768, 18968, 29918, 9638,
                                                                       9728, 20318, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 32018, 0, 3,
                                                                       29918, 19068, 30068, 9728,
                                                                       9818, 20468, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 32243, 0, 3,
                                                                       30218, 19268, 30368, 9998,
                                                                       10088, 20618, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 32468, 0, 3,
                                                                       30368, 19368, 30518,
                                                                       10088, 10178, 20768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 32693, 0, 3,
                                                                       30518, 19468, 30668,
                                                                       10178, 10268, 20918,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 32918, 0, 3,
                                                                       30668, 19568, 30818,
                                                                       10268, 10358, 21068,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 33143, 0, 3,
                                                                       30818, 19668, 30968,
                                                                       10358, 10448, 21218,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 33368, 0, 3,
                                                                       31118, 19868, 31343,
                                                                       10628, 10754, 21368,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 33683, 0, 3,
                                                                       31343, 20018, 31568,
                                                                       10754, 10880, 21578,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 33998, 0, 3,
                                                                       31568, 20168, 31793,
                                                                       10880, 11006, 21788,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 34313, 0, 3,
                                                                       31793, 20318, 32018,
                                                                       11006, 11132, 21998,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 34628, 0, 3,
                                                                       32243, 20618, 32468,
                                                                       11384, 11510, 22208,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 34943, 0, 3,
                                                                       32468, 20768, 32693,
                                                                       11510, 11636, 22418,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 35258, 0, 3,
                                                                       32693, 20918, 32918,
                                                                       11636, 11762, 22628,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 35573, 0, 3,
                                                                       32918, 21068, 33143,
                                                                       11762, 11888, 22838,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 35888, 0, 3,
                                                                       33368, 21368, 33683,
                                                                       12140, 12308, 23048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 36308, 0, 3,
                                                                       33683, 21578, 33998,
                                                                       12308, 12476, 23328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 36728, 0, 3,
                                                                       33998, 21788, 34313,
                                                                       12476, 12644, 23608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 37148, 0, 3,
                                                                       34628, 22208, 34943,
                                                                       12980, 13148, 23888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 37568, 0, 3,
                                                                       34943, 22418, 35258,
                                                                       13148, 13316, 24168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 37988, 0, 3,
                                                                       35258, 22628, 35573,
                                                                       13316, 13484, 24448,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 38408, 0, 3,
                                                                       35888, 23048, 36308,
                                                                       13820, 14036, 24728,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 38948, 0, 3,
                                                                       36308, 23328, 36728,
                                                                       14036, 14252, 25088,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 39488, 0, 3,
                                                                       37148, 23888, 37568,
                                                                       14684, 14900, 25448,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 40028, 0, 3,
                                                                       37568, 24168, 37988,
                                                                       14900, 15116, 25808,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 40568, 0, 3,
                                                                       38408, 24728, 38948,
                                                                       15548, 15818, 26168,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 41243, 0, 3,
                                                                       39488, 25448, 40028,
                                                                       16358, 16628, 26618,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_f_x(buffer, 41918, 28058, 31118, 1, 15, ncols, beta);

                    simdgeo::geom_f_y(buffer, 42068, 28058, 31118, 1, 15, ncols, beta);

                    simdgeo::geom_f_z(buffer, 42218, 28058, 31118, 1, 15, ncols, beta);

                    simdgeo::geom_f_x(buffer, 42368, 28688, 32243, 1, 15, ncols, beta);

                    simdgeo::geom_f_y(buffer, 42518, 28688, 32243, 1, 15, ncols, beta);

                    simdgeo::geom_f_z(buffer, 42668, 28688, 32243, 1, 15, ncols, beta);

                    simdgeo::geom_g_x(buffer, 42818, 29318, 33368, 1, 15, ncols, beta);

                    simdgeo::geom_g_y(buffer, 43043, 29318, 33368, 1, 15, ncols, beta);

                    simdgeo::geom_g_z(buffer, 43268, 29318, 33368, 1, 15, ncols, beta);

                    simdgeo::geom_g_x(buffer, 43493, 30218, 34628, 1, 15, ncols, beta);

                    simdgeo::geom_g_y(buffer, 43718, 30218, 34628, 1, 15, ncols, beta);

                    simdgeo::geom_g_z(buffer, 43943, 30218, 34628, 1, 15, ncols, beta);

                    simdgeo::geom_h_x(buffer, 44168, 31118, 35888, 1, 15, ncols, beta);

                    simdgeo::geom_h_y(buffer, 44483, 31118, 35888, 1, 15, ncols, beta);

                    simdgeo::geom_h_z(buffer, 44798, 31118, 35888, 1, 15, ncols, beta);

                    simdgeo::geom_h_x(buffer, 45113, 32243, 37148, 1, 15, ncols, beta);

                    simdgeo::geom_h_y(buffer, 45428, 32243, 37148, 1, 15, ncols, beta);

                    simdgeo::geom_h_z(buffer, 45743, 32243, 37148, 1, 15, ncols, beta);

                    simdgeo::geom_i_x(buffer, 46058, 33368, 38408, 1, 15, ncols, beta);

                    simdgeo::geom_i_y(buffer, 46478, 33368, 38408, 1, 15, ncols, beta);

                    simdgeo::geom_i_z(buffer, 46898, 33368, 38408, 1, 15, ncols, beta);

                    simdgeo::geom_i_x(buffer, 47318, 34628, 39488, 1, 15, ncols, beta);

                    simdgeo::geom_i_y(buffer, 47738, 34628, 39488, 1, 15, ncols, beta);

                    simdgeo::geom_i_z(buffer, 48158, 34628, 39488, 1, 15, ncols, beta);

                    simdgeo::geom_k_x(buffer, 48578, 35888, 40568, 1, 15, ncols, beta);

                    simdgeo::geom_k_y(buffer, 49118, 35888, 40568, 1, 15, ncols, beta);

                    simdgeo::geom_k_z(buffer, 49658, 35888, 40568, 1, 15, ncols, beta);

                    simdgeo::geom_k_x(buffer, 50198, 37148, 41243, 1, 15, ncols, beta);

                    simdgeo::geom_k_y(buffer, 50738, 37148, 41243, 1, 15, ncols, beta);

                    simdgeo::geom_k_z(buffer, 51278, 37148, 41243, 1, 15, ncols, beta);

                    simdfunc::contract_primitives(buffer, 51818, 41918, 150, ncols);

                    simdfunc::contract_primitives(buffer, 52058, 42068, 150, ncols);

                    simdfunc::contract_primitives(buffer, 52298, 42218, 150, ncols);

                    simdfunc::contract_primitives(buffer, 52538, 29318, 150, ncols);

                    simdfunc::contract_primitives(buffer, 52778, 42368, 150, ncols);

                    simdfunc::contract_primitives(buffer, 53018, 42518, 150, ncols);

                    simdfunc::contract_primitives(buffer, 53258, 42668, 150, ncols);

                    simdfunc::contract_primitives(buffer, 53498, 30218, 150, ncols);

                    simdfunc::contract_primitives(buffer, 53738, 42818, 225, ncols);

                    simdfunc::contract_primitives(buffer, 54098, 43043, 225, ncols);

                    simdfunc::contract_primitives(buffer, 54458, 43268, 225, ncols);

                    simdfunc::contract_primitives(buffer, 54818, 31118, 225, ncols);

                    simdfunc::contract_primitives(buffer, 55178, 43493, 225, ncols);

                    simdfunc::contract_primitives(buffer, 55538, 43718, 225, ncols);

                    simdfunc::contract_primitives(buffer, 55898, 43943, 225, ncols);

                    simdfunc::contract_primitives(buffer, 56258, 32243, 225, ncols);

                    simdfunc::contract_primitives(buffer, 56618, 44168, 315, ncols);

                    simdfunc::contract_primitives(buffer, 57122, 44483, 315, ncols);

                    simdfunc::contract_primitives(buffer, 57626, 44798, 315, ncols);

                    simdfunc::contract_primitives(buffer, 58130, 33368, 315, ncols);

                    simdfunc::contract_primitives(buffer, 58634, 45113, 315, ncols);

                    simdfunc::contract_primitives(buffer, 59138, 45428, 315, ncols);

                    simdfunc::contract_primitives(buffer, 59642, 45743, 315, ncols);

                    simdfunc::contract_primitives(buffer, 60146, 34628, 315, ncols);

                    simdfunc::contract_primitives(buffer, 60650, 46058, 420, ncols);

                    simdfunc::contract_primitives(buffer, 61322, 46478, 420, ncols);

                    simdfunc::contract_primitives(buffer, 61994, 46898, 420, ncols);

                    simdfunc::contract_primitives(buffer, 62666, 35888, 420, ncols);

                    simdfunc::contract_primitives(buffer, 63338, 47318, 420, ncols);

                    simdfunc::contract_primitives(buffer, 64010, 47738, 420, ncols);

                    simdfunc::contract_primitives(buffer, 64682, 48158, 420, ncols);

                    simdfunc::contract_primitives(buffer, 65354, 37148, 420, ncols);

                    simdfunc::contract_primitives(buffer, 66026, 48578, 540, ncols);

                    simdfunc::contract_primitives(buffer, 66890, 49118, 540, ncols);

                    simdfunc::contract_primitives(buffer, 67754, 49658, 540, ncols);

                    simdfunc::contract_primitives(buffer, 68618, 50198, 540, ncols);

                    simdfunc::contract_primitives(buffer, 69482, 50738, 540, ncols);

                    simdfunc::contract_primitives(buffer, 70346, 51278, 540, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 51968, 51818, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 52208, 52058, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 52448, 52298, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 52688, 52538, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 52928, 52778, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 53168, 53018, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 53408, 53258, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 53648, 53498, 10, 1, nmax);

        simdtrf::transform_g_inner(buffer, 53963, 53738, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 54323, 54098, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 54683, 54458, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 55043, 54818, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 55403, 55178, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 55763, 55538, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 56123, 55898, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 56483, 56258, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 56933, 56618, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 57437, 57122, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 57941, 57626, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 58445, 58130, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 58949, 58634, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 59453, 59138, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 59957, 59642, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 60461, 60146, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 61070, 60650, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 61742, 61322, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 62414, 61994, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 63086, 62666, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 63758, 63338, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 64430, 64010, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 65102, 64682, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 65774, 65354, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 66566, 66026, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 67430, 66890, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 68294, 67754, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 69158, 68618, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 70022, 69482, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 70886, 70346, 36, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 71210, 51968, 52688, 53963, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 71480, 52208, 52688, 54323, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 71750, 52448, 52688, 54683, 9,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 72020, 52688, 55043, 9, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 72290, 52928, 53648, 55403, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 72560, 53168, 53648, 55763, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 72830, 53408, 53648, 56123, 9,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 73100, 53648, 56483, 9, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 73370, 53963, 55043, 56933, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 73775, 54323, 55043, 57437, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 74180, 54683, 55043, 57941, 9,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 74585, 55043, 58445, 9, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 74990, 55403, 56483, 58949, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 75395, 55763, 56483, 59453, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 75800, 56123, 56483, 59957, 9,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 76205, 56483, 60461, 9, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 76610, 56933, 58445, 61070, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 77177, 57437, 58445, 61742, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 77744, 57941, 58445, 62414, 9,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 78311, 58445, 63086, 9, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 78878, 58949, 60461, 63758, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 79445, 59453, 60461, 64430, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 80012, 59957, 60461, 65102, 9,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 80579, 60461, 65774, 9, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 81146, 61070, 63086, 66566, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 81902, 61742, 63086, 67430, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 82658, 62414, 63086, 68294, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 83414, 63758, 65774, 69158, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 84170, 64430, 65774, 70022, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 84926, 65102, 65774, 70886, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 85682, 71210, 72020, 73370, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 86222, 71480, 72020, 73775, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 86762, 71750, 72020, 74180, 9,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 87302, 72020, 74585, 9, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 87842, 72290, 73100, 74990, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 88382, 72560, 73100, 75395, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 88922, 72830, 73100, 75800, 9,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 89462, 73100, 76205, 9, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 90002, 73370, 74585, 76610, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 90812, 73775, 74585, 77177, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 91622, 74180, 74585, 77744, 9,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 92432, 74585, 78311, 9, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 93242, 74990, 76205, 78878, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 94052, 75395, 76205, 79445, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 94862, 75800, 76205, 80012, 9,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 95672, 76205, 80579, 9, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 96482, 76610, 78311, 81146, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 97616, 77177, 78311, 81902, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 98750, 77744, 78311, 82658, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 99884, 78878, 80579, 83414, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 101018, 79445, 80579, 84170, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 102152, 80012, 80579, 84926, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 103286, 85682, 87302, 90002, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 104186, 86222, 87302, 90812, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 105086, 86762, 87302, 91622, 9,
                                          nmax);

        simdtrf::compute_hrr_ff(buffer, coordinates, 105986, 87302, 92432, 9, nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 106886, 87842, 89462, 93242, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 107786, 88382, 89462, 94052, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 108686, 88922, 89462, 94862, 9,
                                          nmax);

        simdtrf::compute_hrr_ff(buffer, coordinates, 109586, 89462, 95672, 9, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 110486, 90002, 92432, 96482, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 111836, 90812, 92432, 97616, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 113186, 91622, 92432, 98750, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 114536, 93242, 95672, 99884, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 115886, 94052, 95672, 101018, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 117236, 94862, 95672, 102152, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gf_out_of_second(buffer, coordinates, 118586, 103286,
                                                        105986, 110486, 9, nmax);

        simdtrf::compute_hrr_geom_010y_gf_out_of_second(buffer, coordinates, 119936, 104186,
                                                        105986, 111836, 9, nmax);

        simdtrf::compute_hrr_geom_010z_gf_out_of_second(buffer, coordinates, 121286, 105086,
                                                        105986, 113186, 9, nmax);

        simdtrf::compute_hrr_geom_010x_gf_out_of_second(buffer, coordinates, 122636, 106886,
                                                        109586, 114536, 9, nmax);

        simdtrf::compute_hrr_geom_010y_gf_out_of_second(buffer, coordinates, 123986, 107786,
                                                        109586, 115886, 9, nmax);

        simdtrf::compute_hrr_geom_010z_gf_out_of_second(buffer, coordinates, 125336, 108686,
                                                        109586, 117236, 9, nmax);

        simdtrf::transform_f_inner(buffer, 126686, 122636, 15, 9, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 126686, 63, nmax);

        simdtrf::transform_f_inner(buffer, 126686, 123986, 15, 9, nmax);

        simdtrf::transform_g_outer(values + 567 * nvalues + n * npairs, nvalues, buffer, 126686,
                                   63, nmax);

        simdtrf::transform_f_inner(buffer, 126686, 125336, 15, 9, nmax);

        simdtrf::transform_g_outer(values + 1134 * nvalues + n * npairs, nvalues, buffer, 126686,
                                   63, nmax);

        simdtrf::transform_f_inner(buffer, 126686, 118586, 15, 9, nmax);

        simdtrf::transform_g_outer(values + 1701 * nvalues + n * npairs, nvalues, buffer, 126686,
                                   63, nmax);

        simdtrf::transform_f_inner(buffer, 126686, 119936, 15, 9, nmax);

        simdtrf::transform_g_outer(values + 2268 * nvalues + n * npairs, nvalues, buffer, 126686,
                                   63, nmax);

        simdtrf::transform_f_inner(buffer, 126686, 121286, 15, 9, nmax);

        simdtrf::transform_g_outer(values + 2835 * nvalues + n * npairs, nvalues, buffer, 126686,
                                   63, nmax);
    }

    for (size_t m = 0; m < 3402; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
