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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecDFI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSII.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferGeom010XDF.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010YDF.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010ZDF.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_dfi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_dfi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 83694, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2730 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 83694, 58472, 13093, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1436, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1439, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1442, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1445, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1448, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1451, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1454, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1457, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1460, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1463, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1466, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1469, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1472, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1475, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1478, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1481, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1484, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1487, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1490, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1493, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1496, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1499, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1502, 3, 9, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1511, 3, 10, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1520, 3, 11, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1529, 3, 12, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1538, 3, 13, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1547, 3, 14, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1556, 3, 15, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1565, 3, 16, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1574, 3, 17, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1583, 3, 18, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1592, 3, 23, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1601, 3, 24, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1610, 3, 25, 82,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1619, 3, 26, 85,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1628, 3, 27, 88,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1637, 3, 28, 91,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1646, 3, 29, 94,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1655, 3, 30, 97,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1664, 3, 31, 100,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1673, 3, 32, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1682, 3, 40, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1700, 3, 43, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1718, 3, 46, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1736, 3, 49, 136,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1754, 3, 52, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1772, 3, 55, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1790, 3, 58, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1808, 3, 61, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1826, 3, 64, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1844, 3, 76, 184,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1862, 3, 79, 190,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1880, 3, 82, 196,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1898, 3, 85, 202,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1916, 3, 88, 208,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1934, 3, 91, 214,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1952, 3, 94, 220,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1970, 3, 97, 226,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1988, 3, 100, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2006, 3, 118, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2036, 3, 124, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2066, 3, 130, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2096, 3, 136, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2126, 3, 142, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2156, 3, 148, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2186, 3, 154, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2216, 3, 160, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2246, 3, 184, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2276, 3, 190, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2306, 3, 196, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2336, 3, 202, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2366, 3, 208, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2396, 3, 214, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2426, 3, 220, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2456, 3, 226, 428,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2486, 3, 258, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2531, 3, 268, 483,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2576, 3, 278, 498,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2621, 3, 288, 513,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2666, 3, 298, 528,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2711, 3, 308, 543,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2756, 3, 318, 558,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2801, 3, 358, 603,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2846, 3, 368, 618,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2891, 3, 378, 633,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2936, 3, 388, 648,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2981, 3, 398, 663,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3026, 3, 408, 678,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3071, 3, 418, 693,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3116, 3, 468, 750,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3179, 3, 483, 771,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3242, 3, 498, 792,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3305, 3, 513, 813,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3368, 3, 528, 834,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3431, 3, 543, 855,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3494, 3, 603, 918,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3557, 3, 618, 939,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3620, 3, 633, 960,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3683, 3, 648, 981,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3746, 3, 663,
                                                                       1002, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3809, 3, 678,
                                                                       1023, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3872, 3, 750,
                                                                       1100, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3956, 3, 771,
                                                                       1128, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4040, 3, 792,
                                                                       1156, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4124, 3, 813,
                                                                       1184, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4208, 3, 834,
                                                                       1212, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4292, 3, 918,
                                                                       1296, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4376, 3, 939,
                                                                       1324, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4460, 3, 960,
                                                                       1352, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4544, 3, 981,
                                                                       1380, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4628, 3, 1002,
                                                                       1408, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4712, 3, 7, 8,
                                                                       1436, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4718, 3, 8, 9,
                                                                       1439, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4724, 3, 9, 10,
                                                                       1442, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4730, 3, 10, 11,
                                                                       1445, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4736, 3, 11, 12,
                                                                       1448, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4742, 3, 12, 13,
                                                                       1451, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4748, 3, 13, 14,
                                                                       1454, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4754, 3, 14, 15,
                                                                       1457, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4760, 3, 15, 16,
                                                                       1460, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4766, 3, 16, 17,
                                                                       1463, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4772, 3, 17, 18,
                                                                       1466, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4778, 3, 21, 22,
                                                                       1469, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4784, 3, 22, 23,
                                                                       1472, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4790, 3, 23, 24,
                                                                       1475, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4796, 3, 24, 25,
                                                                       1478, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4802, 3, 25, 26,
                                                                       1481, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4808, 3, 26, 27,
                                                                       1484, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4814, 3, 27, 28,
                                                                       1487, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4820, 3, 28, 29,
                                                                       1490, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4826, 3, 29, 30,
                                                                       1493, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4832, 3, 30, 31,
                                                                       1496, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4838, 3, 31, 32,
                                                                       1499, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4844, 0, 3, 4712,
                                                                       1436, 4718, 1502, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4862, 0, 3, 4718,
                                                                       1439, 4724, 1511, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4880, 0, 3, 4724,
                                                                       1442, 4730, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4898, 0, 3, 4730,
                                                                       1445, 4736, 1529, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4916, 0, 3, 4736,
                                                                       1448, 4742, 1538, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4934, 0, 3, 4742,
                                                                       1451, 4748, 1547, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4952, 0, 3, 4748,
                                                                       1454, 4754, 1556, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4970, 0, 3, 4754,
                                                                       1457, 4760, 1565, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4988, 0, 3, 4760,
                                                                       1460, 4766, 1574, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5006, 0, 3, 4766,
                                                                       1463, 4772, 1583, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5024, 0, 3, 4778,
                                                                       1469, 4784, 1592, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5042, 0, 3, 4784,
                                                                       1472, 4790, 1601, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5060, 0, 3, 4790,
                                                                       1475, 4796, 1610, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5078, 0, 3, 4796,
                                                                       1478, 4802, 1619, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5096, 0, 3, 4802,
                                                                       1481, 4808, 1628, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5114, 0, 3, 4808,
                                                                       1484, 4814, 1637, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5132, 0, 3, 4814,
                                                                       1487, 4820, 1646, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5150, 0, 3, 4820,
                                                                       1490, 4826, 1655, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5168, 0, 3, 4826,
                                                                       1493, 4832, 1664, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5186, 0, 3, 4832,
                                                                       1496, 4838, 1673, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5204, 0, 3, 4844,
                                                                       1502, 4862, 106, 112,
                                                                       1682, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5240, 0, 3, 4862,
                                                                       1511, 4880, 112, 118,
                                                                       1700, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5276, 0, 3, 4880,
                                                                       1520, 4898, 118, 124,
                                                                       1718, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5312, 0, 3, 4898,
                                                                       1529, 4916, 124, 130,
                                                                       1736, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5348, 0, 3, 4916,
                                                                       1538, 4934, 130, 136,
                                                                       1754, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5384, 0, 3, 4934,
                                                                       1547, 4952, 136, 142,
                                                                       1772, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5420, 0, 3, 4952,
                                                                       1556, 4970, 142, 148,
                                                                       1790, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5456, 0, 3, 4970,
                                                                       1565, 4988, 148, 154,
                                                                       1808, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5492, 0, 3, 4988,
                                                                       1574, 5006, 154, 160,
                                                                       1826, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5528, 0, 3, 5024,
                                                                       1592, 5042, 172, 178,
                                                                       1844, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5564, 0, 3, 5042,
                                                                       1601, 5060, 178, 184,
                                                                       1862, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5600, 0, 3, 5060,
                                                                       1610, 5078, 184, 190,
                                                                       1880, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5636, 0, 3, 5078,
                                                                       1619, 5096, 190, 196,
                                                                       1898, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5672, 0, 3, 5096,
                                                                       1628, 5114, 196, 202,
                                                                       1916, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5708, 0, 3, 5114,
                                                                       1637, 5132, 202, 208,
                                                                       1934, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5744, 0, 3, 5132,
                                                                       1646, 5150, 208, 214,
                                                                       1952, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5780, 0, 3, 5150,
                                                                       1655, 5168, 214, 220,
                                                                       1970, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5816, 0, 3, 5168,
                                                                       1664, 5186, 220, 226,
                                                                       1988, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5852, 0, 3, 5204,
                                                                       1682, 5240, 238, 248,
                                                                       2006, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5912, 0, 3, 5240,
                                                                       1700, 5276, 248, 258,
                                                                       2036, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5972, 0, 3, 5276,
                                                                       1718, 5312, 258, 268,
                                                                       2066, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6032, 0, 3, 5312,
                                                                       1736, 5348, 268, 278,
                                                                       2096, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6092, 0, 3, 5348,
                                                                       1754, 5384, 278, 288,
                                                                       2126, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6152, 0, 3, 5384,
                                                                       1772, 5420, 288, 298,
                                                                       2156, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6212, 0, 3, 5420,
                                                                       1790, 5456, 298, 308,
                                                                       2186, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6272, 0, 3, 5456,
                                                                       1808, 5492, 308, 318,
                                                                       2216, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6332, 0, 3, 5528,
                                                                       1844, 5564, 338, 348,
                                                                       2246, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6392, 0, 3, 5564,
                                                                       1862, 5600, 348, 358,
                                                                       2276, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6452, 0, 3, 5600,
                                                                       1880, 5636, 358, 368,
                                                                       2306, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6512, 0, 3, 5636,
                                                                       1898, 5672, 368, 378,
                                                                       2336, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6572, 0, 3, 5672,
                                                                       1916, 5708, 378, 388,
                                                                       2366, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6632, 0, 3, 5708,
                                                                       1934, 5744, 388, 398,
                                                                       2396, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6692, 0, 3, 5744,
                                                                       1952, 5780, 398, 408,
                                                                       2426, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6752, 0, 3, 5780,
                                                                       1970, 5816, 408, 418,
                                                                       2456, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6812, 0, 3, 5852,
                                                                       2006, 5912, 438, 453,
                                                                       2486, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6902, 0, 3, 5912,
                                                                       2036, 5972, 453, 468,
                                                                       2531, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6992, 0, 3, 5972,
                                                                       2066, 6032, 468, 483,
                                                                       2576, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7082, 0, 3, 6032,
                                                                       2096, 6092, 483, 498,
                                                                       2621, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7172, 0, 3, 6092,
                                                                       2126, 6152, 498, 513,
                                                                       2666, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7262, 0, 3, 6152,
                                                                       2156, 6212, 513, 528,
                                                                       2711, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7352, 0, 3, 6212,
                                                                       2186, 6272, 528, 543,
                                                                       2756, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7442, 0, 3, 6332,
                                                                       2246, 6392, 573, 588,
                                                                       2801, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7532, 0, 3, 6392,
                                                                       2276, 6452, 588, 603,
                                                                       2846, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7622, 0, 3, 6452,
                                                                       2306, 6512, 603, 618,
                                                                       2891, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7712, 0, 3, 6512,
                                                                       2336, 6572, 618, 633,
                                                                       2936, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7802, 0, 3, 6572,
                                                                       2366, 6632, 633, 648,
                                                                       2981, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7892, 0, 3, 6632,
                                                                       2396, 6692, 648, 663,
                                                                       3026, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7982, 0, 3, 6692,
                                                                       2426, 6752, 663, 678,
                                                                       3071, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8072, 0, 3, 6812,
                                                                       2486, 6902, 708, 729,
                                                                       3116, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8198, 0, 3, 6902,
                                                                       2531, 6992, 729, 750,
                                                                       3179, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8324, 0, 3, 6992,
                                                                       2576, 7082, 750, 771,
                                                                       3242, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8450, 0, 3, 7082,
                                                                       2621, 7172, 771, 792,
                                                                       3305, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8576, 0, 3, 7172,
                                                                       2666, 7262, 792, 813,
                                                                       3368, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8702, 0, 3, 7262,
                                                                       2711, 7352, 813, 834,
                                                                       3431, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8828, 0, 3, 7442,
                                                                       2801, 7532, 876, 897,
                                                                       3494, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8954, 0, 3, 7532,
                                                                       2846, 7622, 897, 918,
                                                                       3557, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9080, 0, 3, 7622,
                                                                       2891, 7712, 918, 939,
                                                                       3620, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9206, 0, 3, 7712,
                                                                       2936, 7802, 939, 960,
                                                                       3683, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9332, 0, 3, 7802,
                                                                       2981, 7892, 960, 981,
                                                                       3746, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9458, 0, 3, 7892,
                                                                       3026, 7982, 981, 1002,
                                                                       3809, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9584, 0, 3, 8072,
                                                                       3116, 8198, 1044, 1072,
                                                                       3872, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9752, 0, 3, 8198,
                                                                       3179, 8324, 1072, 1100,
                                                                       3956, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9920, 0, 3, 8324,
                                                                       3242, 8450, 1100, 1128,
                                                                       4040, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10088, 0, 3, 8450,
                                                                       3305, 8576, 1128, 1156,
                                                                       4124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10256, 0, 3, 8576,
                                                                       3368, 8702, 1156, 1184,
                                                                       4208, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10424, 0, 3, 8828,
                                                                       3494, 8954, 1240, 1268,
                                                                       4292, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10592, 0, 3, 8954,
                                                                       3557, 9080, 1268, 1296,
                                                                       4376, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10760, 0, 3, 9080,
                                                                       3620, 9206, 1296, 1324,
                                                                       4460, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10928, 0, 3, 9206,
                                                                       3683, 9332, 1324, 1352,
                                                                       4544, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11096, 0, 3, 9332,
                                                                       3746, 9458, 1352, 1380,
                                                                       4628, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11264, 3, 1436,
                                                                       1439, 4724, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11274, 3, 1439,
                                                                       1442, 4730, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11284, 3, 1442,
                                                                       1445, 4736, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11294, 3, 1445,
                                                                       1448, 4742, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11304, 3, 1448,
                                                                       1451, 4748, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11314, 3, 1451,
                                                                       1454, 4754, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11324, 3, 1454,
                                                                       1457, 4760, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11334, 3, 1457,
                                                                       1460, 4766, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11344, 3, 1460,
                                                                       1463, 4772, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11354, 3, 1469,
                                                                       1472, 4790, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11364, 3, 1472,
                                                                       1475, 4796, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11374, 3, 1475,
                                                                       1478, 4802, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11384, 3, 1478,
                                                                       1481, 4808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11394, 3, 1481,
                                                                       1484, 4814, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11404, 3, 1484,
                                                                       1487, 4820, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11414, 3, 1487,
                                                                       1490, 4826, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11424, 3, 1490,
                                                                       1493, 4832, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11434, 3, 1493,
                                                                       1496, 4838, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11444, 0, 3,
                                                                       11264, 4724, 11274, 4880,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11474, 0, 3,
                                                                       11274, 4730, 11284, 4898,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11504, 0, 3,
                                                                       11284, 4736, 11294, 4916,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11534, 0, 3,
                                                                       11294, 4742, 11304, 4934,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11564, 0, 3,
                                                                       11304, 4748, 11314, 4952,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11594, 0, 3,
                                                                       11314, 4754, 11324, 4970,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11624, 0, 3,
                                                                       11324, 4760, 11334, 4988,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11654, 0, 3,
                                                                       11334, 4766, 11344, 5006,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11684, 0, 3,
                                                                       11354, 4790, 11364, 5060,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11714, 0, 3,
                                                                       11364, 4796, 11374, 5078,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11744, 0, 3,
                                                                       11374, 4802, 11384, 5096,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11774, 0, 3,
                                                                       11384, 4808, 11394, 5114,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11804, 0, 3,
                                                                       11394, 4814, 11404, 5132,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11834, 0, 3,
                                                                       11404, 4820, 11414, 5150,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11864, 0, 3,
                                                                       11414, 4826, 11424, 5168,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11894, 0, 3,
                                                                       11424, 4832, 11434, 5186,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11924, 0, 3,
                                                                       11444, 4880, 11474, 1682,
                                                                       1700, 5276, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11984, 0, 3,
                                                                       11474, 4898, 11504, 1700,
                                                                       1718, 5312, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12044, 0, 3,
                                                                       11504, 4916, 11534, 1718,
                                                                       1736, 5348, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12104, 0, 3,
                                                                       11534, 4934, 11564, 1736,
                                                                       1754, 5384, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12164, 0, 3,
                                                                       11564, 4952, 11594, 1754,
                                                                       1772, 5420, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12224, 0, 3,
                                                                       11594, 4970, 11624, 1772,
                                                                       1790, 5456, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12284, 0, 3,
                                                                       11624, 4988, 11654, 1790,
                                                                       1808, 5492, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12344, 0, 3,
                                                                       11684, 5060, 11714, 1844,
                                                                       1862, 5600, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12404, 0, 3,
                                                                       11714, 5078, 11744, 1862,
                                                                       1880, 5636, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12464, 0, 3,
                                                                       11744, 5096, 11774, 1880,
                                                                       1898, 5672, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12524, 0, 3,
                                                                       11774, 5114, 11804, 1898,
                                                                       1916, 5708, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12584, 0, 3,
                                                                       11804, 5132, 11834, 1916,
                                                                       1934, 5744, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12644, 0, 3,
                                                                       11834, 5150, 11864, 1934,
                                                                       1952, 5780, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12704, 0, 3,
                                                                       11864, 5168, 11894, 1952,
                                                                       1970, 5816, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12764, 0, 3,
                                                                       11924, 5276, 11984, 2006,
                                                                       2036, 5972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12864, 0, 3,
                                                                       11984, 5312, 12044, 2036,
                                                                       2066, 6032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12964, 0, 3,
                                                                       12044, 5348, 12104, 2066,
                                                                       2096, 6092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13064, 0, 3,
                                                                       12104, 5384, 12164, 2096,
                                                                       2126, 6152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13164, 0, 3,
                                                                       12164, 5420, 12224, 2126,
                                                                       2156, 6212, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13264, 0, 3,
                                                                       12224, 5456, 12284, 2156,
                                                                       2186, 6272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13364, 0, 3,
                                                                       12344, 5600, 12404, 2246,
                                                                       2276, 6452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13464, 0, 3,
                                                                       12404, 5636, 12464, 2276,
                                                                       2306, 6512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13564, 0, 3,
                                                                       12464, 5672, 12524, 2306,
                                                                       2336, 6572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13664, 0, 3,
                                                                       12524, 5708, 12584, 2336,
                                                                       2366, 6632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13764, 0, 3,
                                                                       12584, 5744, 12644, 2366,
                                                                       2396, 6692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13864, 0, 3,
                                                                       12644, 5780, 12704, 2396,
                                                                       2426, 6752, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13964, 0, 3,
                                                                       12764, 5972, 12864, 2486,
                                                                       2531, 6992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14114, 0, 3,
                                                                       12864, 6032, 12964, 2531,
                                                                       2576, 7082, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14264, 0, 3,
                                                                       12964, 6092, 13064, 2576,
                                                                       2621, 7172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14414, 0, 3,
                                                                       13064, 6152, 13164, 2621,
                                                                       2666, 7262, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14564, 0, 3,
                                                                       13164, 6212, 13264, 2666,
                                                                       2711, 7352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14714, 0, 3,
                                                                       13364, 6452, 13464, 2801,
                                                                       2846, 7622, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14864, 0, 3,
                                                                       13464, 6512, 13564, 2846,
                                                                       2891, 7712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15014, 0, 3,
                                                                       13564, 6572, 13664, 2891,
                                                                       2936, 7802, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15164, 0, 3,
                                                                       13664, 6632, 13764, 2936,
                                                                       2981, 7892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15314, 0, 3,
                                                                       13764, 6692, 13864, 2981,
                                                                       3026, 7982, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15464, 0, 3,
                                                                       13964, 6992, 14114, 3116,
                                                                       3179, 8324, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15674, 0, 3,
                                                                       14114, 7082, 14264, 3179,
                                                                       3242, 8450, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15884, 0, 3,
                                                                       14264, 7172, 14414, 3242,
                                                                       3305, 8576, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16094, 0, 3,
                                                                       14414, 7262, 14564, 3305,
                                                                       3368, 8702, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16304, 0, 3,
                                                                       14714, 7622, 14864, 3494,
                                                                       3557, 9080, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16514, 0, 3,
                                                                       14864, 7712, 15014, 3557,
                                                                       3620, 9206, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16724, 0, 3,
                                                                       15014, 7802, 15164, 3620,
                                                                       3683, 9332, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16934, 0, 3,
                                                                       15164, 7892, 15314, 3683,
                                                                       3746, 9458, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 17144, 0, 3,
                                                                       15464, 8324, 15674, 3872,
                                                                       3956, 9920, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 17424, 0, 3,
                                                                       15674, 8450, 15884, 3956,
                                                                       4040, 10088, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 17704, 0, 3,
                                                                       15884, 8576, 16094, 4040,
                                                                       4124, 10256, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 17984, 0, 3,
                                                                       16304, 9080, 16514, 4292,
                                                                       4376, 10760, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18264, 0, 3,
                                                                       16514, 9206, 16724, 4376,
                                                                       4460, 10928, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 18544, 0, 3,
                                                                       16724, 9332, 16934, 4460,
                                                                       4544, 11096, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18824, 3, 4712,
                                                                       4718, 11264, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18839, 3, 4718,
                                                                       4724, 11274, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18854, 3, 4724,
                                                                       4730, 11284, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18869, 3, 4730,
                                                                       4736, 11294, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18884, 3, 4736,
                                                                       4742, 11304, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18899, 3, 4742,
                                                                       4748, 11314, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18914, 3, 4748,
                                                                       4754, 11324, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18929, 3, 4754,
                                                                       4760, 11334, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18944, 3, 4760,
                                                                       4766, 11344, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18959, 3, 4778,
                                                                       4784, 11354, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18974, 3, 4784,
                                                                       4790, 11364, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18989, 3, 4790,
                                                                       4796, 11374, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 19004, 3, 4796,
                                                                       4802, 11384, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 19019, 3, 4802,
                                                                       4808, 11394, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 19034, 3, 4808,
                                                                       4814, 11404, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 19049, 3, 4814,
                                                                       4820, 11414, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 19064, 3, 4820,
                                                                       4826, 11424, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 19079, 3, 4826,
                                                                       4832, 11434, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19094, 0, 3,
                                                                       18824, 11264, 18839, 4844,
                                                                       4862, 11444, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19139, 0, 3,
                                                                       18839, 11274, 18854, 4862,
                                                                       4880, 11474, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19184, 0, 3,
                                                                       18854, 11284, 18869, 4880,
                                                                       4898, 11504, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19229, 0, 3,
                                                                       18869, 11294, 18884, 4898,
                                                                       4916, 11534, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19274, 0, 3,
                                                                       18884, 11304, 18899, 4916,
                                                                       4934, 11564, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19319, 0, 3,
                                                                       18899, 11314, 18914, 4934,
                                                                       4952, 11594, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19364, 0, 3,
                                                                       18914, 11324, 18929, 4952,
                                                                       4970, 11624, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19409, 0, 3,
                                                                       18929, 11334, 18944, 4970,
                                                                       4988, 11654, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19454, 0, 3,
                                                                       18959, 11354, 18974, 5024,
                                                                       5042, 11684, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19499, 0, 3,
                                                                       18974, 11364, 18989, 5042,
                                                                       5060, 11714, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19544, 0, 3,
                                                                       18989, 11374, 19004, 5060,
                                                                       5078, 11744, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19589, 0, 3,
                                                                       19004, 11384, 19019, 5078,
                                                                       5096, 11774, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19634, 0, 3,
                                                                       19019, 11394, 19034, 5096,
                                                                       5114, 11804, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19679, 0, 3,
                                                                       19034, 11404, 19049, 5114,
                                                                       5132, 11834, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19724, 0, 3,
                                                                       19049, 11414, 19064, 5132,
                                                                       5150, 11864, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19769, 0, 3,
                                                                       19064, 11424, 19079, 5150,
                                                                       5168, 11894, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19814, 0, 3,
                                                                       19094, 11444, 19139, 5204,
                                                                       5240, 11924, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19904, 0, 3,
                                                                       19139, 11474, 19184, 5240,
                                                                       5276, 11984, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19994, 0, 3,
                                                                       19184, 11504, 19229, 5276,
                                                                       5312, 12044, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20084, 0, 3,
                                                                       19229, 11534, 19274, 5312,
                                                                       5348, 12104, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20174, 0, 3,
                                                                       19274, 11564, 19319, 5348,
                                                                       5384, 12164, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20264, 0, 3,
                                                                       19319, 11594, 19364, 5384,
                                                                       5420, 12224, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20354, 0, 3,
                                                                       19364, 11624, 19409, 5420,
                                                                       5456, 12284, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20444, 0, 3,
                                                                       19454, 11684, 19499, 5528,
                                                                       5564, 12344, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20534, 0, 3,
                                                                       19499, 11714, 19544, 5564,
                                                                       5600, 12404, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20624, 0, 3,
                                                                       19544, 11744, 19589, 5600,
                                                                       5636, 12464, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20714, 0, 3,
                                                                       19589, 11774, 19634, 5636,
                                                                       5672, 12524, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20804, 0, 3,
                                                                       19634, 11804, 19679, 5672,
                                                                       5708, 12584, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20894, 0, 3,
                                                                       19679, 11834, 19724, 5708,
                                                                       5744, 12644, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20984, 0, 3,
                                                                       19724, 11864, 19769, 5744,
                                                                       5780, 12704, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21074, 0, 3,
                                                                       19814, 11924, 19904, 5852,
                                                                       5912, 12764, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21224, 0, 3,
                                                                       19904, 11984, 19994, 5912,
                                                                       5972, 12864, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21374, 0, 3,
                                                                       19994, 12044, 20084, 5972,
                                                                       6032, 12964, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21524, 0, 3,
                                                                       20084, 12104, 20174, 6032,
                                                                       6092, 13064, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21674, 0, 3,
                                                                       20174, 12164, 20264, 6092,
                                                                       6152, 13164, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21824, 0, 3,
                                                                       20264, 12224, 20354, 6152,
                                                                       6212, 13264, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21974, 0, 3,
                                                                       20444, 12344, 20534, 6332,
                                                                       6392, 13364, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 22124, 0, 3,
                                                                       20534, 12404, 20624, 6392,
                                                                       6452, 13464, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 22274, 0, 3,
                                                                       20624, 12464, 20714, 6452,
                                                                       6512, 13564, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 22424, 0, 3,
                                                                       20714, 12524, 20804, 6512,
                                                                       6572, 13664, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 22574, 0, 3,
                                                                       20804, 12584, 20894, 6572,
                                                                       6632, 13764, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 22724, 0, 3,
                                                                       20894, 12644, 20984, 6632,
                                                                       6692, 13864, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22874, 0, 3,
                                                                       21074, 12764, 21224, 6812,
                                                                       6902, 13964, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23099, 0, 3,
                                                                       21224, 12864, 21374, 6902,
                                                                       6992, 14114, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23324, 0, 3,
                                                                       21374, 12964, 21524, 6992,
                                                                       7082, 14264, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23549, 0, 3,
                                                                       21524, 13064, 21674, 7082,
                                                                       7172, 14414, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23774, 0, 3,
                                                                       21674, 13164, 21824, 7172,
                                                                       7262, 14564, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23999, 0, 3,
                                                                       21974, 13364, 22124, 7442,
                                                                       7532, 14714, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 24224, 0, 3,
                                                                       22124, 13464, 22274, 7532,
                                                                       7622, 14864, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 24449, 0, 3,
                                                                       22274, 13564, 22424, 7622,
                                                                       7712, 15014, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 24674, 0, 3,
                                                                       22424, 13664, 22574, 7712,
                                                                       7802, 15164, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 24899, 0, 3,
                                                                       22574, 13764, 22724, 7802,
                                                                       7892, 15314, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 25124, 0, 3,
                                                                       22874, 13964, 23099, 8072,
                                                                       8198, 15464, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 25439, 0, 3,
                                                                       23099, 14114, 23324, 8198,
                                                                       8324, 15674, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 25754, 0, 3,
                                                                       23324, 14264, 23549, 8324,
                                                                       8450, 15884, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 26069, 0, 3,
                                                                       23549, 14414, 23774, 8450,
                                                                       8576, 16094, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 26384, 0, 3,
                                                                       23999, 14714, 24224, 8828,
                                                                       8954, 16304, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 26699, 0, 3,
                                                                       24224, 14864, 24449, 8954,
                                                                       9080, 16514, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 27014, 0, 3,
                                                                       24449, 15014, 24674, 9080,
                                                                       9206, 16724, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 27329, 0, 3,
                                                                       24674, 15164, 24899, 9206,
                                                                       9332, 16934, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 27644, 0, 3,
                                                                       25124, 15464, 25439, 9584,
                                                                       9752, 17144, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 28064, 0, 3,
                                                                       25439, 15674, 25754, 9752,
                                                                       9920, 17424, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 28484, 0, 3,
                                                                       25754, 15884, 26069, 9920,
                                                                       10088, 17704, ncols,
                                                                       gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 28904, 0, 3,
                                                                       26384, 16304, 26699,
                                                                       10424, 10592, 17984,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 29324, 0, 3,
                                                                       26699, 16514, 27014,
                                                                       10592, 10760, 18264,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 29744, 0, 3,
                                                                       27014, 16724, 27329,
                                                                       10760, 10928, 18544,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30164, 3, 11264,
                                                                       11274, 18854, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30185, 3, 11274,
                                                                       11284, 18869, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30206, 3, 11284,
                                                                       11294, 18884, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30227, 3, 11294,
                                                                       11304, 18899, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30248, 3, 11304,
                                                                       11314, 18914, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30269, 3, 11314,
                                                                       11324, 18929, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30290, 3, 11324,
                                                                       11334, 18944, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30311, 3, 11354,
                                                                       11364, 18989, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30332, 3, 11364,
                                                                       11374, 19004, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30353, 3, 11374,
                                                                       11384, 19019, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30374, 3, 11384,
                                                                       11394, 19034, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30395, 3, 11394,
                                                                       11404, 19049, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30416, 3, 11404,
                                                                       11414, 19064, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30437, 3, 11414,
                                                                       11424, 19079, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 30458, 0, 3,
                                                                       30164, 18854, 30185,
                                                                       11444, 11474, 19184,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 30521, 0, 3,
                                                                       30185, 18869, 30206,
                                                                       11474, 11504, 19229,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 30584, 0, 3,
                                                                       30206, 18884, 30227,
                                                                       11504, 11534, 19274,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 30647, 0, 3,
                                                                       30227, 18899, 30248,
                                                                       11534, 11564, 19319,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 30710, 0, 3,
                                                                       30248, 18914, 30269,
                                                                       11564, 11594, 19364,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 30773, 0, 3,
                                                                       30269, 18929, 30290,
                                                                       11594, 11624, 19409,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 30836, 0, 3,
                                                                       30311, 18989, 30332,
                                                                       11684, 11714, 19544,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 30899, 0, 3,
                                                                       30332, 19004, 30353,
                                                                       11714, 11744, 19589,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 30962, 0, 3,
                                                                       30353, 19019, 30374,
                                                                       11744, 11774, 19634,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31025, 0, 3,
                                                                       30374, 19034, 30395,
                                                                       11774, 11804, 19679,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31088, 0, 3,
                                                                       30395, 19049, 30416,
                                                                       11804, 11834, 19724,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31151, 0, 3,
                                                                       30416, 19064, 30437,
                                                                       11834, 11864, 19769,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 31214, 0, 3,
                                                                       30458, 19184, 30521,
                                                                       11924, 11984, 19994,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 31340, 0, 3,
                                                                       30521, 19229, 30584,
                                                                       11984, 12044, 20084,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 31466, 0, 3,
                                                                       30584, 19274, 30647,
                                                                       12044, 12104, 20174,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 31592, 0, 3,
                                                                       30647, 19319, 30710,
                                                                       12104, 12164, 20264,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 31718, 0, 3,
                                                                       30710, 19364, 30773,
                                                                       12164, 12224, 20354,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 31844, 0, 3,
                                                                       30836, 19544, 30899,
                                                                       12344, 12404, 20624,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 31970, 0, 3,
                                                                       30899, 19589, 30962,
                                                                       12404, 12464, 20714,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32096, 0, 3,
                                                                       30962, 19634, 31025,
                                                                       12464, 12524, 20804,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32222, 0, 3,
                                                                       31025, 19679, 31088,
                                                                       12524, 12584, 20894,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32348, 0, 3,
                                                                       31088, 19724, 31151,
                                                                       12584, 12644, 20984,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 32474, 0, 3,
                                                                       31214, 19994, 31340,
                                                                       12764, 12864, 21374,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 32684, 0, 3,
                                                                       31340, 20084, 31466,
                                                                       12864, 12964, 21524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 32894, 0, 3,
                                                                       31466, 20174, 31592,
                                                                       12964, 13064, 21674,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 33104, 0, 3,
                                                                       31592, 20264, 31718,
                                                                       13064, 13164, 21824,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 33314, 0, 3,
                                                                       31844, 20624, 31970,
                                                                       13364, 13464, 22274,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 33524, 0, 3,
                                                                       31970, 20714, 32096,
                                                                       13464, 13564, 22424,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 33734, 0, 3,
                                                                       32096, 20804, 32222,
                                                                       13564, 13664, 22574,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 33944, 0, 3,
                                                                       32222, 20894, 32348,
                                                                       13664, 13764, 22724,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 34154, 0, 3,
                                                                       32474, 21374, 32684,
                                                                       13964, 14114, 23324,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 34469, 0, 3,
                                                                       32684, 21524, 32894,
                                                                       14114, 14264, 23549,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 34784, 0, 3,
                                                                       32894, 21674, 33104,
                                                                       14264, 14414, 23774,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 35099, 0, 3,
                                                                       33314, 22274, 33524,
                                                                       14714, 14864, 24449,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 35414, 0, 3,
                                                                       33524, 22424, 33734,
                                                                       14864, 15014, 24674,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 35729, 0, 3,
                                                                       33734, 22574, 33944,
                                                                       15014, 15164, 24899,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 36044, 0, 3,
                                                                       34154, 23324, 34469,
                                                                       15464, 15674, 25754,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 36485, 0, 3,
                                                                       34469, 23549, 34784,
                                                                       15674, 15884, 26069,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 36926, 0, 3,
                                                                       35099, 24449, 35414,
                                                                       16304, 16514, 27014,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 37367, 0, 3,
                                                                       35414, 24674, 35729,
                                                                       16514, 16724, 27329,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 37808, 0, 3,
                                                                       36044, 25754, 36485,
                                                                       17144, 17424, 28484,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 38396, 0, 3,
                                                                       36926, 27014, 37367,
                                                                       17984, 18264, 29744,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 38984, 3, 18824,
                                                                       18839, 30164, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 39012, 3, 18839,
                                                                       18854, 30185, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 39040, 3, 18854,
                                                                       18869, 30206, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 39068, 3, 18869,
                                                                       18884, 30227, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 39096, 3, 18884,
                                                                       18899, 30248, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 39124, 3, 18899,
                                                                       18914, 30269, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 39152, 3, 18914,
                                                                       18929, 30290, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 39180, 3, 18959,
                                                                       18974, 30311, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 39208, 3, 18974,
                                                                       18989, 30332, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 39236, 3, 18989,
                                                                       19004, 30353, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 39264, 3, 19004,
                                                                       19019, 30374, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 39292, 3, 19019,
                                                                       19034, 30395, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 39320, 3, 19034,
                                                                       19049, 30416, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 39348, 3, 19049,
                                                                       19064, 30437, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 39376, 0, 3,
                                                                       38984, 30164, 39012,
                                                                       19094, 19139, 30458,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 39460, 0, 3,
                                                                       39012, 30185, 39040,
                                                                       19139, 19184, 30521,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 39544, 0, 3,
                                                                       39040, 30206, 39068,
                                                                       19184, 19229, 30584,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 39628, 0, 3,
                                                                       39068, 30227, 39096,
                                                                       19229, 19274, 30647,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 39712, 0, 3,
                                                                       39096, 30248, 39124,
                                                                       19274, 19319, 30710,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 39796, 0, 3,
                                                                       39124, 30269, 39152,
                                                                       19319, 19364, 30773,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 39880, 0, 3,
                                                                       39180, 30311, 39208,
                                                                       19454, 19499, 30836,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 39964, 0, 3,
                                                                       39208, 30332, 39236,
                                                                       19499, 19544, 30899,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 40048, 0, 3,
                                                                       39236, 30353, 39264,
                                                                       19544, 19589, 30962,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 40132, 0, 3,
                                                                       39264, 30374, 39292,
                                                                       19589, 19634, 31025,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 40216, 0, 3,
                                                                       39292, 30395, 39320,
                                                                       19634, 19679, 31088,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 40300, 0, 3,
                                                                       39320, 30416, 39348,
                                                                       19679, 19724, 31151,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 40384, 0, 3,
                                                                       39376, 30458, 39460,
                                                                       19814, 19904, 31214,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 40552, 0, 3,
                                                                       39460, 30521, 39544,
                                                                       19904, 19994, 31340,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 40720, 0, 3,
                                                                       39544, 30584, 39628,
                                                                       19994, 20084, 31466,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 40888, 0, 3,
                                                                       39628, 30647, 39712,
                                                                       20084, 20174, 31592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 41056, 0, 3,
                                                                       39712, 30710, 39796,
                                                                       20174, 20264, 31718,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 41224, 0, 3,
                                                                       39880, 30836, 39964,
                                                                       20444, 20534, 31844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 41392, 0, 3,
                                                                       39964, 30899, 40048,
                                                                       20534, 20624, 31970,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 41560, 0, 3,
                                                                       40048, 30962, 40132,
                                                                       20624, 20714, 32096,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 41728, 0, 3,
                                                                       40132, 31025, 40216,
                                                                       20714, 20804, 32222,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 41896, 0, 3,
                                                                       40216, 31088, 40300,
                                                                       20804, 20894, 32348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 42064, 0, 3,
                                                                       40384, 31214, 40552,
                                                                       21074, 21224, 32474,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 42344, 0, 3,
                                                                       40552, 31340, 40720,
                                                                       21224, 21374, 32684,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 42624, 0, 3,
                                                                       40720, 31466, 40888,
                                                                       21374, 21524, 32894,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 42904, 0, 3,
                                                                       40888, 31592, 41056,
                                                                       21524, 21674, 33104,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 43184, 0, 3,
                                                                       41224, 31844, 41392,
                                                                       21974, 22124, 33314,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 43464, 0, 3,
                                                                       41392, 31970, 41560,
                                                                       22124, 22274, 33524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 43744, 0, 3,
                                                                       41560, 32096, 41728,
                                                                       22274, 22424, 33734,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 44024, 0, 3,
                                                                       41728, 32222, 41896,
                                                                       22424, 22574, 33944,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 44304, 0, 3,
                                                                       42064, 32474, 42344,
                                                                       22874, 23099, 34154,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 44724, 0, 3,
                                                                       42344, 32684, 42624,
                                                                       23099, 23324, 34469,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 45144, 0, 3,
                                                                       42624, 32894, 42904,
                                                                       23324, 23549, 34784,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 45564, 0, 3,
                                                                       43184, 33314, 43464,
                                                                       23999, 24224, 35099,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 45984, 0, 3,
                                                                       43464, 33524, 43744,
                                                                       24224, 24449, 35414,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 46404, 0, 3,
                                                                       43744, 33734, 44024,
                                                                       24449, 24674, 35729,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 46824, 0, 3,
                                                                       44304, 34154, 44724,
                                                                       25124, 25439, 36044,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 47412, 0, 3,
                                                                       44724, 34469, 45144,
                                                                       25439, 25754, 36485,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 48000, 0, 3,
                                                                       45564, 35099, 45984,
                                                                       26384, 26699, 36926,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 48588, 0, 3,
                                                                       45984, 35414, 46404,
                                                                       26699, 27014, 37367,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 49176, 0, 3,
                                                                       46824, 36044, 47412,
                                                                       27644, 28064, 37808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 49960, 0, 3,
                                                                       48000, 36926, 48588,
                                                                       28904, 29324, 38396,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_f_x(buffer, 50744, 40384, 44304, 1, 28, ncols, beta);

                    simdgeo::geom_f_y(buffer, 51024, 40384, 44304, 1, 28, ncols, beta);

                    simdgeo::geom_f_z(buffer, 51304, 40384, 44304, 1, 28, ncols, beta);

                    simdgeo::geom_f_x(buffer, 51584, 41224, 45564, 1, 28, ncols, beta);

                    simdgeo::geom_f_y(buffer, 51864, 41224, 45564, 1, 28, ncols, beta);

                    simdgeo::geom_f_z(buffer, 52144, 41224, 45564, 1, 28, ncols, beta);

                    simdgeo::geom_g_x(buffer, 52424, 42064, 46824, 1, 28, ncols, beta);

                    simdgeo::geom_g_y(buffer, 52844, 42064, 46824, 1, 28, ncols, beta);

                    simdgeo::geom_g_z(buffer, 53264, 42064, 46824, 1, 28, ncols, beta);

                    simdgeo::geom_g_x(buffer, 53684, 43184, 48000, 1, 28, ncols, beta);

                    simdgeo::geom_g_y(buffer, 54104, 43184, 48000, 1, 28, ncols, beta);

                    simdgeo::geom_g_z(buffer, 54524, 43184, 48000, 1, 28, ncols, beta);

                    simdgeo::geom_h_x(buffer, 54944, 44304, 49176, 1, 28, ncols, beta);

                    simdgeo::geom_h_y(buffer, 55532, 44304, 49176, 1, 28, ncols, beta);

                    simdgeo::geom_h_z(buffer, 56120, 44304, 49176, 1, 28, ncols, beta);

                    simdgeo::geom_h_x(buffer, 56708, 45564, 49960, 1, 28, ncols, beta);

                    simdgeo::geom_h_y(buffer, 57296, 45564, 49960, 1, 28, ncols, beta);

                    simdgeo::geom_h_z(buffer, 57884, 45564, 49960, 1, 28, ncols, beta);

                    simdfunc::contract_primitives(buffer, 58472, 50744, 280, ncols);

                    simdfunc::contract_primitives(buffer, 58882, 51024, 280, ncols);

                    simdfunc::contract_primitives(buffer, 59292, 51304, 280, ncols);

                    simdfunc::contract_primitives(buffer, 59702, 42064, 280, ncols);

                    simdfunc::contract_primitives(buffer, 60112, 51584, 280, ncols);

                    simdfunc::contract_primitives(buffer, 60522, 51864, 280, ncols);

                    simdfunc::contract_primitives(buffer, 60932, 52144, 280, ncols);

                    simdfunc::contract_primitives(buffer, 61342, 43184, 280, ncols);

                    simdfunc::contract_primitives(buffer, 61752, 52424, 420, ncols);

                    simdfunc::contract_primitives(buffer, 62367, 52844, 420, ncols);

                    simdfunc::contract_primitives(buffer, 62982, 53264, 420, ncols);

                    simdfunc::contract_primitives(buffer, 63597, 44304, 420, ncols);

                    simdfunc::contract_primitives(buffer, 64212, 53684, 420, ncols);

                    simdfunc::contract_primitives(buffer, 64827, 54104, 420, ncols);

                    simdfunc::contract_primitives(buffer, 65442, 54524, 420, ncols);

                    simdfunc::contract_primitives(buffer, 66057, 45564, 420, ncols);

                    simdfunc::contract_primitives(buffer, 66672, 54944, 588, ncols);

                    simdfunc::contract_primitives(buffer, 67533, 55532, 588, ncols);

                    simdfunc::contract_primitives(buffer, 68394, 56120, 588, ncols);

                    simdfunc::contract_primitives(buffer, 69255, 56708, 588, ncols);

                    simdfunc::contract_primitives(buffer, 70116, 57296, 588, ncols);

                    simdfunc::contract_primitives(buffer, 70977, 57884, 588, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 58752, 58472, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 59162, 58882, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 59572, 59292, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 59982, 59702, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 60392, 60112, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 60802, 60522, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 61212, 60932, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 61622, 61342, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 62172, 61752, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 62787, 62367, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 63402, 62982, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 64017, 63597, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 64632, 64212, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 65247, 64827, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 65862, 65442, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 66477, 66057, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 67260, 66672, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 68121, 67533, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 68982, 68394, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 69843, 69255, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 70704, 70116, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 71565, 70977, 21, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 71838, 58752, 59982, 62172, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 72228, 59162, 59982, 62787, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 72618, 59572, 59982, 63402, 13,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 73008, 59982, 64017, 13, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 73398, 60392, 61622, 64632, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 73788, 60802, 61622, 65247, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 74178, 61212, 61622, 65862, 13,
                                          nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 74568, 61622, 66477, 13, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 74958, 62172, 64017, 67260, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 75543, 62787, 64017, 68121, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 76128, 63402, 64017, 68982, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 76713, 64632, 66477, 69843, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 77298, 65247, 66477, 70704, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 77883, 65862, 66477, 71565, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 78468, 71838, 73008, 74958, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 79248, 72228, 73008, 75543, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 80028, 72618, 73008, 76128, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 80808, 73398, 74568, 76713, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 81588, 73788, 74568, 77298, 13,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 82368, 74178, 74568, 77883, 13,
                                          nmax);

        simdtrf::transform_f_inner(buffer, 83148, 80808, 6, 13, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 83148, 91, nmax);

        simdtrf::transform_f_inner(buffer, 83148, 81588, 6, 13, nmax);

        simdtrf::transform_d_outer(values + 455 * nvalues + n * npairs, nvalues, buffer, 83148,
                                   91, nmax);

        simdtrf::transform_f_inner(buffer, 83148, 82368, 6, 13, nmax);

        simdtrf::transform_d_outer(values + 910 * nvalues + n * npairs, nvalues, buffer, 83148,
                                   91, nmax);

        simdtrf::transform_f_inner(buffer, 83148, 78468, 6, 13, nmax);

        simdtrf::transform_d_outer(values + 1365 * nvalues + n * npairs, nvalues, buffer, 83148,
                                   91, nmax);

        simdtrf::transform_f_inner(buffer, 83148, 79248, 6, 13, nmax);

        simdtrf::transform_d_outer(values + 1820 * nvalues + n * npairs, nvalues, buffer, 83148,
                                   91, nmax);

        simdtrf::transform_f_inner(buffer, 83148, 80028, 6, 13, nmax);

        simdtrf::transform_d_outer(values + 2275 * nvalues + n * npairs, nvalues, buffer, 83148,
                                   91, nmax);
    }

    for (size_t m = 0; m < 2730; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
