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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecDGH.hpp"

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
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_dgh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_dgh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 85760, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2970 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 85760, 56516, 14284, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10, 11, 12}, ncols,
                                                            fj, i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 19, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 56, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 59, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 62, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 65, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 68, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 71, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 74, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 77, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 80, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 83, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 86, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 89, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 92, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 95, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 98, 0, 3, 7, 8,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 104, 0, 3, 8, 9,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 110, 0, 3, 9, 10,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 116, 0, 3, 10, 11,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 122, 0, 3, 11, 12,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 128, 0, 3, 12, 13,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 134, 0, 3, 13, 14,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 140, 0, 3, 14, 15,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 146, 0, 3, 15, 16,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 152, 0, 3, 16, 17,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 158, 0, 3, 20, 21,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 164, 0, 3, 21, 22,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 170, 0, 3, 22, 23,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 176, 0, 3, 23, 24,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 182, 0, 3, 24, 25,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 188, 0, 3, 25, 26,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 194, 0, 3, 26, 27,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 200, 0, 3, 27, 28,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 206, 0, 3, 28, 29,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 212, 0, 3, 29, 30,
                                                                       92, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 218, 0, 3, 32, 35,
                                                                       98, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 35, 38,
                                                                       104, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 38, 41,
                                                                       110, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 248, 0, 3, 41, 44,
                                                                       116, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 44, 47,
                                                                       122, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 47, 50,
                                                                       128, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 50, 53,
                                                                       134, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 53, 56,
                                                                       140, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 56, 59,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 65, 68,
                                                                       158, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 68, 71,
                                                                       164, 170, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 71, 74,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 74, 77,
                                                                       176, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 77, 80,
                                                                       182, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 80, 83,
                                                                       188, 194, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 83, 86,
                                                                       194, 200, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 86, 89,
                                                                       200, 206, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 89, 92,
                                                                       206, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 98,
                                                                       104, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 413, 0, 3, 104,
                                                                       110, 228, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 110,
                                                                       116, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 443, 0, 3, 116,
                                                                       122, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 122,
                                                                       128, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 473, 0, 3, 128,
                                                                       134, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 488, 0, 3, 134,
                                                                       140, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 503, 0, 3, 140,
                                                                       146, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 158,
                                                                       164, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 533, 0, 3, 164,
                                                                       170, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 548, 0, 3, 170,
                                                                       176, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 563, 0, 3, 176,
                                                                       182, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 578, 0, 3, 182,
                                                                       188, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 593, 0, 3, 188,
                                                                       194, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 608, 0, 3, 194,
                                                                       200, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 623, 0, 3, 200,
                                                                       206, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 638, 0, 3, 218,
                                                                       228, 398, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 659, 0, 3, 228,
                                                                       238, 413, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 680, 0, 3, 238,
                                                                       248, 428, 443, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 701, 0, 3, 248,
                                                                       258, 443, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 722, 0, 3, 258,
                                                                       268, 458, 473, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 743, 0, 3, 268,
                                                                       278, 473, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 764, 0, 3, 278,
                                                                       288, 488, 503, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 785, 0, 3, 308,
                                                                       318, 518, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 806, 0, 3, 318,
                                                                       328, 533, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 827, 0, 3, 328,
                                                                       338, 548, 563, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 848, 0, 3, 338,
                                                                       348, 563, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 869, 0, 3, 348,
                                                                       358, 578, 593, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 890, 0, 3, 358,
                                                                       368, 593, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 911, 0, 3, 368,
                                                                       378, 608, 623, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 932, 0, 3, 398,
                                                                       413, 638, 659, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 960, 0, 3, 413,
                                                                       428, 659, 680, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 988, 0, 3, 428,
                                                                       443, 680, 701, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1016, 0, 3, 443,
                                                                       458, 701, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 458,
                                                                       473, 722, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 473,
                                                                       488, 743, 764, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 518,
                                                                       533, 785, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 533,
                                                                       548, 806, 827, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 548,
                                                                       563, 827, 848, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 563,
                                                                       578, 848, 869, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 578,
                                                                       593, 869, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 593,
                                                                       608, 890, 911, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 638,
                                                                       659, 932, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1304, 0, 3, 659,
                                                                       680, 960, 988, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1340, 0, 3, 680,
                                                                       701, 988, 1016, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1376, 0, 3, 701,
                                                                       722, 1016, 1044, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1412, 0, 3, 722,
                                                                       743, 1044, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1448, 0, 3, 785,
                                                                       806, 1100, 1128, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1484, 0, 3, 806,
                                                                       827, 1128, 1156, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 827,
                                                                       848, 1156, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1556, 0, 3, 848,
                                                                       869, 1184, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1592, 0, 3, 869,
                                                                       890, 1212, 1240, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1628, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1631, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1634, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1637, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1640, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1643, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1646, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1649, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1652, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1655, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1658, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1661, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1664, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1667, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1670, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1673, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1676, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1679, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1682, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1685, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1688, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1691, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1694, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1697, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1700, 3, 9, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1709, 3, 10, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1718, 3, 11, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1727, 3, 12, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1736, 3, 13, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1745, 3, 14, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1754, 3, 15, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1763, 3, 16, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1772, 3, 17, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1781, 3, 22, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1790, 3, 23, 74,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1799, 3, 24, 77,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1808, 3, 25, 80,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1817, 3, 26, 83,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1826, 3, 27, 86,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1835, 3, 28, 89,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1844, 3, 29, 92,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1853, 3, 30, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1862, 3, 32, 98,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1880, 3, 35, 104,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1898, 3, 38, 110,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1916, 3, 41, 116,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1934, 3, 44, 122,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1952, 3, 47, 128,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1970, 3, 50, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1988, 3, 53, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2006, 3, 56, 146,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2024, 3, 59, 152,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2042, 3, 65, 158,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2060, 3, 68, 164,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2078, 3, 71, 170,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2096, 3, 74, 176,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2114, 3, 77, 182,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2132, 3, 80, 188,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2150, 3, 83, 194,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2168, 3, 86, 200,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2186, 3, 89, 206,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2204, 3, 92, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2222, 3, 98, 218,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2252, 3, 104, 228,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2282, 3, 110, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2312, 3, 116, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2342, 3, 122, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2372, 3, 128, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2402, 3, 134, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2432, 3, 140, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2462, 3, 146, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2492, 3, 158, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2522, 3, 164, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2552, 3, 170, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2582, 3, 176, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2612, 3, 182, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2642, 3, 188, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2672, 3, 194, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2702, 3, 200, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2732, 3, 206, 388,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2762, 3, 218, 398,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2807, 3, 228, 413,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2852, 3, 238, 428,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2897, 3, 248, 443,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2942, 3, 258, 458,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2987, 3, 268, 473,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3032, 3, 278, 488,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3077, 3, 288, 503,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3122, 3, 308, 518,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3167, 3, 318, 533,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3212, 3, 328, 548,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3257, 3, 338, 563,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3302, 3, 348, 578,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3347, 3, 358, 593,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3392, 3, 368, 608,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3437, 3, 378, 623,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3482, 3, 398, 638,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3545, 3, 413, 659,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3608, 3, 428, 680,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3671, 3, 443, 701,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3734, 3, 458, 722,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3797, 3, 473, 743,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3860, 3, 488, 764,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3923, 3, 518, 785,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3986, 3, 533, 806,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4049, 3, 548, 827,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4112, 3, 563, 848,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4175, 3, 578, 869,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4238, 3, 593, 890,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4301, 3, 608, 911,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4364, 3, 638, 932,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4448, 3, 659, 960,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4532, 3, 680, 988,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4616, 3, 701,
                                                                       1016, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4700, 3, 722,
                                                                       1044, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4784, 3, 743,
                                                                       1072, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4868, 3, 785,
                                                                       1100, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4952, 3, 806,
                                                                       1128, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5036, 3, 827,
                                                                       1156, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5120, 3, 848,
                                                                       1184, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5204, 3, 869,
                                                                       1212, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5288, 3, 890,
                                                                       1240, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5372, 3, 932,
                                                                       1268, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5480, 3, 960,
                                                                       1304, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5588, 3, 988,
                                                                       1340, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5696, 3, 1016,
                                                                       1376, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5804, 3, 1044,
                                                                       1412, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5912, 3, 1100,
                                                                       1448, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6020, 3, 1128,
                                                                       1484, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6128, 3, 1156,
                                                                       1520, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6236, 3, 1184,
                                                                       1556, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6344, 3, 1212,
                                                                       1592, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6452, 3, 7, 8,
                                                                       1634, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6458, 3, 8, 9,
                                                                       1637, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6464, 3, 9, 10,
                                                                       1640, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6470, 3, 10, 11,
                                                                       1643, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6476, 3, 11, 12,
                                                                       1646, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6482, 3, 12, 13,
                                                                       1649, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6488, 3, 13, 14,
                                                                       1652, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6494, 3, 14, 15,
                                                                       1655, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6500, 3, 15, 16,
                                                                       1658, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6506, 3, 16, 17,
                                                                       1661, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6512, 3, 20, 21,
                                                                       1670, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6518, 3, 21, 22,
                                                                       1673, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6524, 3, 22, 23,
                                                                       1676, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6530, 3, 23, 24,
                                                                       1679, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6536, 3, 24, 25,
                                                                       1682, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6542, 3, 25, 26,
                                                                       1685, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6548, 3, 26, 27,
                                                                       1688, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6554, 3, 27, 28,
                                                                       1691, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6560, 3, 28, 29,
                                                                       1694, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6566, 3, 29, 30,
                                                                       1697, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6572, 0, 3, 6452,
                                                                       1634, 6458, 1700, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6590, 0, 3, 6458,
                                                                       1637, 6464, 1709, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6608, 0, 3, 6464,
                                                                       1640, 6470, 1718, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6626, 0, 3, 6470,
                                                                       1643, 6476, 1727, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6644, 0, 3, 6476,
                                                                       1646, 6482, 1736, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6662, 0, 3, 6482,
                                                                       1649, 6488, 1745, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6680, 0, 3, 6488,
                                                                       1652, 6494, 1754, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6698, 0, 3, 6494,
                                                                       1655, 6500, 1763, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6716, 0, 3, 6500,
                                                                       1658, 6506, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6734, 0, 3, 6512,
                                                                       1670, 6518, 1781, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6752, 0, 3, 6518,
                                                                       1673, 6524, 1790, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6770, 0, 3, 6524,
                                                                       1676, 6530, 1799, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6788, 0, 3, 6530,
                                                                       1679, 6536, 1808, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6806, 0, 3, 6536,
                                                                       1682, 6542, 1817, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6824, 0, 3, 6542,
                                                                       1685, 6548, 1826, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6842, 0, 3, 6548,
                                                                       1688, 6554, 1835, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6860, 0, 3, 6554,
                                                                       1691, 6560, 1844, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6878, 0, 3, 6560,
                                                                       1694, 6566, 1853, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6896, 0, 3, 6572,
                                                                       1700, 6590, 98, 104, 1898,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6932, 0, 3, 6590,
                                                                       1709, 6608, 104, 110,
                                                                       1916, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6968, 0, 3, 6608,
                                                                       1718, 6626, 110, 116,
                                                                       1934, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7004, 0, 3, 6626,
                                                                       1727, 6644, 116, 122,
                                                                       1952, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7040, 0, 3, 6644,
                                                                       1736, 6662, 122, 128,
                                                                       1970, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7076, 0, 3, 6662,
                                                                       1745, 6680, 128, 134,
                                                                       1988, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7112, 0, 3, 6680,
                                                                       1754, 6698, 134, 140,
                                                                       2006, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7148, 0, 3, 6698,
                                                                       1763, 6716, 140, 146,
                                                                       2024, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7184, 0, 3, 6734,
                                                                       1781, 6752, 158, 164,
                                                                       2078, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7220, 0, 3, 6752,
                                                                       1790, 6770, 164, 170,
                                                                       2096, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7256, 0, 3, 6770,
                                                                       1799, 6788, 170, 176,
                                                                       2114, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7292, 0, 3, 6788,
                                                                       1808, 6806, 176, 182,
                                                                       2132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7328, 0, 3, 6806,
                                                                       1817, 6824, 182, 188,
                                                                       2150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7364, 0, 3, 6824,
                                                                       1826, 6842, 188, 194,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 6842,
                                                                       1835, 6860, 194, 200,
                                                                       2186, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7436, 0, 3, 6860,
                                                                       1844, 6878, 200, 206,
                                                                       2204, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7472, 0, 3, 6896,
                                                                       1898, 6932, 218, 228,
                                                                       2282, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7532, 0, 3, 6932,
                                                                       1916, 6968, 228, 238,
                                                                       2312, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7592, 0, 3, 6968,
                                                                       1934, 7004, 238, 248,
                                                                       2342, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7652, 0, 3, 7004,
                                                                       1952, 7040, 248, 258,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7712, 0, 3, 7040,
                                                                       1970, 7076, 258, 268,
                                                                       2402, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7772, 0, 3, 7076,
                                                                       1988, 7112, 268, 278,
                                                                       2432, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7832, 0, 3, 7112,
                                                                       2006, 7148, 278, 288,
                                                                       2462, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7892, 0, 3, 7184,
                                                                       2078, 7220, 308, 318,
                                                                       2552, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7952, 0, 3, 7220,
                                                                       2096, 7256, 318, 328,
                                                                       2582, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8012, 0, 3, 7256,
                                                                       2114, 7292, 328, 338,
                                                                       2612, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8072, 0, 3, 7292,
                                                                       2132, 7328, 338, 348,
                                                                       2642, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8132, 0, 3, 7328,
                                                                       2150, 7364, 348, 358,
                                                                       2672, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8192, 0, 3, 7364,
                                                                       2168, 7400, 358, 368,
                                                                       2702, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8252, 0, 3, 7400,
                                                                       2186, 7436, 368, 378,
                                                                       2732, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8312, 0, 3, 7472,
                                                                       2282, 7532, 398, 413,
                                                                       2852, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8402, 0, 3, 7532,
                                                                       2312, 7592, 413, 428,
                                                                       2897, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8492, 0, 3, 7592,
                                                                       2342, 7652, 428, 443,
                                                                       2942, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8582, 0, 3, 7652,
                                                                       2372, 7712, 443, 458,
                                                                       2987, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8672, 0, 3, 7712,
                                                                       2402, 7772, 458, 473,
                                                                       3032, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8762, 0, 3, 7772,
                                                                       2432, 7832, 473, 488,
                                                                       3077, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8852, 0, 3, 7892,
                                                                       2552, 7952, 518, 533,
                                                                       3212, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8942, 0, 3, 7952,
                                                                       2582, 8012, 533, 548,
                                                                       3257, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9032, 0, 3, 8012,
                                                                       2612, 8072, 548, 563,
                                                                       3302, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9122, 0, 3, 8072,
                                                                       2642, 8132, 563, 578,
                                                                       3347, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9212, 0, 3, 8132,
                                                                       2672, 8192, 578, 593,
                                                                       3392, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9302, 0, 3, 8192,
                                                                       2702, 8252, 593, 608,
                                                                       3437, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9392, 0, 3, 8312,
                                                                       2852, 8402, 638, 659,
                                                                       3608, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9518, 0, 3, 8402,
                                                                       2897, 8492, 659, 680,
                                                                       3671, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9644, 0, 3, 8492,
                                                                       2942, 8582, 680, 701,
                                                                       3734, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9770, 0, 3, 8582,
                                                                       2987, 8672, 701, 722,
                                                                       3797, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9896, 0, 3, 8672,
                                                                       3032, 8762, 722, 743,
                                                                       3860, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10022, 0, 3, 8852,
                                                                       3212, 8942, 785, 806,
                                                                       4049, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10148, 0, 3, 8942,
                                                                       3257, 9032, 806, 827,
                                                                       4112, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10274, 0, 3, 9032,
                                                                       3302, 9122, 827, 848,
                                                                       4175, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10400, 0, 3, 9122,
                                                                       3347, 9212, 848, 869,
                                                                       4238, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10526, 0, 3, 9212,
                                                                       3392, 9302, 869, 890,
                                                                       4301, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10652, 0, 3, 9392,
                                                                       3608, 9518, 932, 960,
                                                                       4532, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10820, 0, 3, 9518,
                                                                       3671, 9644, 960, 988,
                                                                       4616, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10988, 0, 3, 9644,
                                                                       3734, 9770, 988, 1016,
                                                                       4700, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11156, 0, 3, 9770,
                                                                       3797, 9896, 1016, 1044,
                                                                       4784, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11324, 0, 3,
                                                                       10022, 4049, 10148, 1100,
                                                                       1128, 5036, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11492, 0, 3,
                                                                       10148, 4112, 10274, 1128,
                                                                       1156, 5120, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11660, 0, 3,
                                                                       10274, 4175, 10400, 1156,
                                                                       1184, 5204, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11828, 0, 3,
                                                                       10400, 4238, 10526, 1184,
                                                                       1212, 5288, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11996, 0, 3,
                                                                       10652, 4532, 10820, 1268,
                                                                       1304, 5588, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12212, 0, 3,
                                                                       10820, 4616, 10988, 1304,
                                                                       1340, 5696, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12428, 0, 3,
                                                                       10988, 4700, 11156, 1340,
                                                                       1376, 5804, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12644, 0, 3,
                                                                       11324, 5036, 11492, 1448,
                                                                       1484, 6128, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12860, 0, 3,
                                                                       11492, 5120, 11660, 1484,
                                                                       1520, 6236, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13076, 0, 3,
                                                                       11660, 5204, 11828, 1520,
                                                                       1556, 6344, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13292, 3, 1628,
                                                                       1631, 6452, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13302, 3, 1631,
                                                                       1634, 6458, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13312, 3, 1634,
                                                                       1637, 6464, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13322, 3, 1637,
                                                                       1640, 6470, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13332, 3, 1640,
                                                                       1643, 6476, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13342, 3, 1643,
                                                                       1646, 6482, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13352, 3, 1646,
                                                                       1649, 6488, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13362, 3, 1649,
                                                                       1652, 6494, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13372, 3, 1652,
                                                                       1655, 6500, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13382, 3, 1655,
                                                                       1658, 6506, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13392, 3, 1664,
                                                                       1667, 6512, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13402, 3, 1667,
                                                                       1670, 6518, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13412, 3, 1670,
                                                                       1673, 6524, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13422, 3, 1673,
                                                                       1676, 6530, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13432, 3, 1676,
                                                                       1679, 6536, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13442, 3, 1679,
                                                                       1682, 6542, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13452, 3, 1682,
                                                                       1685, 6548, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13462, 3, 1685,
                                                                       1688, 6554, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13472, 3, 1688,
                                                                       1691, 6560, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 13482, 3, 1691,
                                                                       1694, 6566, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13492, 0, 3,
                                                                       13292, 6452, 13302, 6572,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13522, 0, 3,
                                                                       13302, 6458, 13312, 6590,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13552, 0, 3,
                                                                       13312, 6464, 13322, 6608,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13582, 0, 3,
                                                                       13322, 6470, 13332, 6626,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13612, 0, 3,
                                                                       13332, 6476, 13342, 6644,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13642, 0, 3,
                                                                       13342, 6482, 13352, 6662,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13672, 0, 3,
                                                                       13352, 6488, 13362, 6680,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13702, 0, 3,
                                                                       13362, 6494, 13372, 6698,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13732, 0, 3,
                                                                       13372, 6500, 13382, 6716,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13762, 0, 3,
                                                                       13392, 6512, 13402, 6734,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13792, 0, 3,
                                                                       13402, 6518, 13412, 6752,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13822, 0, 3,
                                                                       13412, 6524, 13422, 6770,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13852, 0, 3,
                                                                       13422, 6530, 13432, 6788,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13882, 0, 3,
                                                                       13432, 6536, 13442, 6806,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13912, 0, 3,
                                                                       13442, 6542, 13452, 6824,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13942, 0, 3,
                                                                       13452, 6548, 13462, 6842,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13972, 0, 3,
                                                                       13462, 6554, 13472, 6860,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14002, 0, 3,
                                                                       13472, 6560, 13482, 6878,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14032, 0, 3,
                                                                       13492, 6572, 13522, 1862,
                                                                       1880, 6896, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14092, 0, 3,
                                                                       13522, 6590, 13552, 1880,
                                                                       1898, 6932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14152, 0, 3,
                                                                       13552, 6608, 13582, 1898,
                                                                       1916, 6968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14212, 0, 3,
                                                                       13582, 6626, 13612, 1916,
                                                                       1934, 7004, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14272, 0, 3,
                                                                       13612, 6644, 13642, 1934,
                                                                       1952, 7040, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14332, 0, 3,
                                                                       13642, 6662, 13672, 1952,
                                                                       1970, 7076, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14392, 0, 3,
                                                                       13672, 6680, 13702, 1970,
                                                                       1988, 7112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14452, 0, 3,
                                                                       13702, 6698, 13732, 1988,
                                                                       2006, 7148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14512, 0, 3,
                                                                       13762, 6734, 13792, 2042,
                                                                       2060, 7184, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14572, 0, 3,
                                                                       13792, 6752, 13822, 2060,
                                                                       2078, 7220, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14632, 0, 3,
                                                                       13822, 6770, 13852, 2078,
                                                                       2096, 7256, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14692, 0, 3,
                                                                       13852, 6788, 13882, 2096,
                                                                       2114, 7292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14752, 0, 3,
                                                                       13882, 6806, 13912, 2114,
                                                                       2132, 7328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14812, 0, 3,
                                                                       13912, 6824, 13942, 2132,
                                                                       2150, 7364, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14872, 0, 3,
                                                                       13942, 6842, 13972, 2150,
                                                                       2168, 7400, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14932, 0, 3,
                                                                       13972, 6860, 14002, 2168,
                                                                       2186, 7436, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14992, 0, 3,
                                                                       14032, 6896, 14092, 2222,
                                                                       2252, 7472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15092, 0, 3,
                                                                       14092, 6932, 14152, 2252,
                                                                       2282, 7532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15192, 0, 3,
                                                                       14152, 6968, 14212, 2282,
                                                                       2312, 7592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15292, 0, 3,
                                                                       14212, 7004, 14272, 2312,
                                                                       2342, 7652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15392, 0, 3,
                                                                       14272, 7040, 14332, 2342,
                                                                       2372, 7712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15492, 0, 3,
                                                                       14332, 7076, 14392, 2372,
                                                                       2402, 7772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15592, 0, 3,
                                                                       14392, 7112, 14452, 2402,
                                                                       2432, 7832, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15692, 0, 3,
                                                                       14512, 7184, 14572, 2492,
                                                                       2522, 7892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15792, 0, 3,
                                                                       14572, 7220, 14632, 2522,
                                                                       2552, 7952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15892, 0, 3,
                                                                       14632, 7256, 14692, 2552,
                                                                       2582, 8012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15992, 0, 3,
                                                                       14692, 7292, 14752, 2582,
                                                                       2612, 8072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 16092, 0, 3,
                                                                       14752, 7328, 14812, 2612,
                                                                       2642, 8132, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 16192, 0, 3,
                                                                       14812, 7364, 14872, 2642,
                                                                       2672, 8192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 16292, 0, 3,
                                                                       14872, 7400, 14932, 2672,
                                                                       2702, 8252, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16392, 0, 3,
                                                                       14992, 7472, 15092, 2762,
                                                                       2807, 8312, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16542, 0, 3,
                                                                       15092, 7532, 15192, 2807,
                                                                       2852, 8402, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16692, 0, 3,
                                                                       15192, 7592, 15292, 2852,
                                                                       2897, 8492, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16842, 0, 3,
                                                                       15292, 7652, 15392, 2897,
                                                                       2942, 8582, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16992, 0, 3,
                                                                       15392, 7712, 15492, 2942,
                                                                       2987, 8672, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17142, 0, 3,
                                                                       15492, 7772, 15592, 2987,
                                                                       3032, 8762, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17292, 0, 3,
                                                                       15692, 7892, 15792, 3122,
                                                                       3167, 8852, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17442, 0, 3,
                                                                       15792, 7952, 15892, 3167,
                                                                       3212, 8942, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17592, 0, 3,
                                                                       15892, 8012, 15992, 3212,
                                                                       3257, 9032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17742, 0, 3,
                                                                       15992, 8072, 16092, 3257,
                                                                       3302, 9122, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17892, 0, 3,
                                                                       16092, 8132, 16192, 3302,
                                                                       3347, 9212, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18042, 0, 3,
                                                                       16192, 8192, 16292, 3347,
                                                                       3392, 9302, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 18192, 0, 3,
                                                                       16392, 8312, 16542, 3482,
                                                                       3545, 9392, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 18402, 0, 3,
                                                                       16542, 8402, 16692, 3545,
                                                                       3608, 9518, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 18612, 0, 3,
                                                                       16692, 8492, 16842, 3608,
                                                                       3671, 9644, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 18822, 0, 3,
                                                                       16842, 8582, 16992, 3671,
                                                                       3734, 9770, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19032, 0, 3,
                                                                       16992, 8672, 17142, 3734,
                                                                       3797, 9896, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19242, 0, 3,
                                                                       17292, 8852, 17442, 3923,
                                                                       3986, 10022, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19452, 0, 3,
                                                                       17442, 8942, 17592, 3986,
                                                                       4049, 10148, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19662, 0, 3,
                                                                       17592, 9032, 17742, 4049,
                                                                       4112, 10274, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19872, 0, 3,
                                                                       17742, 9122, 17892, 4112,
                                                                       4175, 10400, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20082, 0, 3,
                                                                       17892, 9212, 18042, 4175,
                                                                       4238, 10526, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 20292, 0, 3,
                                                                       18192, 9392, 18402, 4364,
                                                                       4448, 10652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 20572, 0, 3,
                                                                       18402, 9518, 18612, 4448,
                                                                       4532, 10820, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 20852, 0, 3,
                                                                       18612, 9644, 18822, 4532,
                                                                       4616, 10988, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21132, 0, 3,
                                                                       18822, 9770, 19032, 4616,
                                                                       4700, 11156, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21412, 0, 3,
                                                                       19242, 10022, 19452, 4868,
                                                                       4952, 11324, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21692, 0, 3,
                                                                       19452, 10148, 19662, 4952,
                                                                       5036, 11492, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21972, 0, 3,
                                                                       19662, 10274, 19872, 5036,
                                                                       5120, 11660, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22252, 0, 3,
                                                                       19872, 10400, 20082, 5120,
                                                                       5204, 11828, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 22532, 0, 3,
                                                                       20292, 10652, 20572, 5372,
                                                                       5480, 11996, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 22892, 0, 3,
                                                                       20572, 10820, 20852, 5480,
                                                                       5588, 12212, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23252, 0, 3,
                                                                       20852, 10988, 21132, 5588,
                                                                       5696, 12428, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23612, 0, 3,
                                                                       21412, 11324, 21692, 5912,
                                                                       6020, 12644, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23972, 0, 3,
                                                                       21692, 11492, 21972, 6020,
                                                                       6128, 12860, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 24332, 0, 3,
                                                                       21972, 11660, 22252, 6128,
                                                                       6236, 13076, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24692, 3, 6452,
                                                                       6458, 13312, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24707, 3, 6458,
                                                                       6464, 13322, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24722, 3, 6464,
                                                                       6470, 13332, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24737, 3, 6470,
                                                                       6476, 13342, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24752, 3, 6476,
                                                                       6482, 13352, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24767, 3, 6482,
                                                                       6488, 13362, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24782, 3, 6488,
                                                                       6494, 13372, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24797, 3, 6494,
                                                                       6500, 13382, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24812, 3, 6512,
                                                                       6518, 13412, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24827, 3, 6518,
                                                                       6524, 13422, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24842, 3, 6524,
                                                                       6530, 13432, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24857, 3, 6530,
                                                                       6536, 13442, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24872, 3, 6536,
                                                                       6542, 13452, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24887, 3, 6542,
                                                                       6548, 13462, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24902, 3, 6548,
                                                                       6554, 13472, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24917, 3, 6554,
                                                                       6560, 13482, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24932, 0, 3,
                                                                       24692, 13312, 24707, 6572,
                                                                       6590, 13552, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24977, 0, 3,
                                                                       24707, 13322, 24722, 6590,
                                                                       6608, 13582, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25022, 0, 3,
                                                                       24722, 13332, 24737, 6608,
                                                                       6626, 13612, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25067, 0, 3,
                                                                       24737, 13342, 24752, 6626,
                                                                       6644, 13642, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25112, 0, 3,
                                                                       24752, 13352, 24767, 6644,
                                                                       6662, 13672, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25157, 0, 3,
                                                                       24767, 13362, 24782, 6662,
                                                                       6680, 13702, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25202, 0, 3,
                                                                       24782, 13372, 24797, 6680,
                                                                       6698, 13732, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25247, 0, 3,
                                                                       24812, 13412, 24827, 6734,
                                                                       6752, 13822, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25292, 0, 3,
                                                                       24827, 13422, 24842, 6752,
                                                                       6770, 13852, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25337, 0, 3,
                                                                       24842, 13432, 24857, 6770,
                                                                       6788, 13882, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25382, 0, 3,
                                                                       24857, 13442, 24872, 6788,
                                                                       6806, 13912, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25427, 0, 3,
                                                                       24872, 13452, 24887, 6806,
                                                                       6824, 13942, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25472, 0, 3,
                                                                       24887, 13462, 24902, 6824,
                                                                       6842, 13972, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 25517, 0, 3,
                                                                       24902, 13472, 24917, 6842,
                                                                       6860, 14002, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25562, 0, 3,
                                                                       24932, 13552, 24977, 6896,
                                                                       6932, 14152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25652, 0, 3,
                                                                       24977, 13582, 25022, 6932,
                                                                       6968, 14212, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25742, 0, 3,
                                                                       25022, 13612, 25067, 6968,
                                                                       7004, 14272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25832, 0, 3,
                                                                       25067, 13642, 25112, 7004,
                                                                       7040, 14332, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25922, 0, 3,
                                                                       25112, 13672, 25157, 7040,
                                                                       7076, 14392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26012, 0, 3,
                                                                       25157, 13702, 25202, 7076,
                                                                       7112, 14452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26102, 0, 3,
                                                                       25247, 13822, 25292, 7184,
                                                                       7220, 14632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26192, 0, 3,
                                                                       25292, 13852, 25337, 7220,
                                                                       7256, 14692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26282, 0, 3,
                                                                       25337, 13882, 25382, 7256,
                                                                       7292, 14752, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26372, 0, 3,
                                                                       25382, 13912, 25427, 7292,
                                                                       7328, 14812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26462, 0, 3,
                                                                       25427, 13942, 25472, 7328,
                                                                       7364, 14872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26552, 0, 3,
                                                                       25472, 13972, 25517, 7364,
                                                                       7400, 14932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26642, 0, 3,
                                                                       25562, 14152, 25652, 7472,
                                                                       7532, 15192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26792, 0, 3,
                                                                       25652, 14212, 25742, 7532,
                                                                       7592, 15292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26942, 0, 3,
                                                                       25742, 14272, 25832, 7592,
                                                                       7652, 15392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27092, 0, 3,
                                                                       25832, 14332, 25922, 7652,
                                                                       7712, 15492, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27242, 0, 3,
                                                                       25922, 14392, 26012, 7712,
                                                                       7772, 15592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27392, 0, 3,
                                                                       26102, 14632, 26192, 7892,
                                                                       7952, 15892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27542, 0, 3,
                                                                       26192, 14692, 26282, 7952,
                                                                       8012, 15992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27692, 0, 3,
                                                                       26282, 14752, 26372, 8012,
                                                                       8072, 16092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27842, 0, 3,
                                                                       26372, 14812, 26462, 8072,
                                                                       8132, 16192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27992, 0, 3,
                                                                       26462, 14872, 26552, 8132,
                                                                       8192, 16292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28142, 0, 3,
                                                                       26642, 15192, 26792, 8312,
                                                                       8402, 16692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28367, 0, 3,
                                                                       26792, 15292, 26942, 8402,
                                                                       8492, 16842, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28592, 0, 3,
                                                                       26942, 15392, 27092, 8492,
                                                                       8582, 16992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28817, 0, 3,
                                                                       27092, 15492, 27242, 8582,
                                                                       8672, 17142, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 29042, 0, 3,
                                                                       27392, 15892, 27542, 8852,
                                                                       8942, 17592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 29267, 0, 3,
                                                                       27542, 15992, 27692, 8942,
                                                                       9032, 17742, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 29492, 0, 3,
                                                                       27692, 16092, 27842, 9032,
                                                                       9122, 17892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 29717, 0, 3,
                                                                       27842, 16192, 27992, 9122,
                                                                       9212, 18042, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 29942, 0, 3,
                                                                       28142, 16692, 28367, 9392,
                                                                       9518, 18612, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 30257, 0, 3,
                                                                       28367, 16842, 28592, 9518,
                                                                       9644, 18822, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 30572, 0, 3,
                                                                       28592, 16992, 28817, 9644,
                                                                       9770, 19032, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 30887, 0, 3,
                                                                       29042, 17592, 29267,
                                                                       10022, 10148, 19662,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 31202, 0, 3,
                                                                       29267, 17742, 29492,
                                                                       10148, 10274, 19872,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 31517, 0, 3,
                                                                       29492, 17892, 29717,
                                                                       10274, 10400, 20082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 31832, 0, 3,
                                                                       29942, 18612, 30257,
                                                                       10652, 10820, 20852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 32252, 0, 3,
                                                                       30257, 18822, 30572,
                                                                       10820, 10988, 21132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 32672, 0, 3,
                                                                       30887, 19662, 31202,
                                                                       11324, 11492, 21972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 33092, 0, 3,
                                                                       31202, 19872, 31517,
                                                                       11492, 11660, 22252,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 33512, 0, 3,
                                                                       31832, 20852, 32252,
                                                                       11996, 12212, 23252,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 34052, 0, 3,
                                                                       32672, 21972, 33092,
                                                                       12644, 12860, 24332,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34592, 3, 13292,
                                                                       13302, 24692, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34613, 3, 13302,
                                                                       13312, 24707, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34634, 3, 13312,
                                                                       13322, 24722, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34655, 3, 13322,
                                                                       13332, 24737, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34676, 3, 13332,
                                                                       13342, 24752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34697, 3, 13342,
                                                                       13352, 24767, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34718, 3, 13352,
                                                                       13362, 24782, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34739, 3, 13362,
                                                                       13372, 24797, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34760, 3, 13392,
                                                                       13402, 24812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34781, 3, 13402,
                                                                       13412, 24827, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34802, 3, 13412,
                                                                       13422, 24842, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34823, 3, 13422,
                                                                       13432, 24857, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34844, 3, 13432,
                                                                       13442, 24872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34865, 3, 13442,
                                                                       13452, 24887, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34886, 3, 13452,
                                                                       13462, 24902, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34907, 3, 13462,
                                                                       13472, 24917, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 34928, 0, 3,
                                                                       34592, 24692, 34613,
                                                                       13492, 13522, 24932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 34991, 0, 3,
                                                                       34613, 24707, 34634,
                                                                       13522, 13552, 24977,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35054, 0, 3,
                                                                       34634, 24722, 34655,
                                                                       13552, 13582, 25022,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35117, 0, 3,
                                                                       34655, 24737, 34676,
                                                                       13582, 13612, 25067,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35180, 0, 3,
                                                                       34676, 24752, 34697,
                                                                       13612, 13642, 25112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35243, 0, 3,
                                                                       34697, 24767, 34718,
                                                                       13642, 13672, 25157,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35306, 0, 3,
                                                                       34718, 24782, 34739,
                                                                       13672, 13702, 25202,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35369, 0, 3,
                                                                       34760, 24812, 34781,
                                                                       13762, 13792, 25247,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35432, 0, 3,
                                                                       34781, 24827, 34802,
                                                                       13792, 13822, 25292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35495, 0, 3,
                                                                       34802, 24842, 34823,
                                                                       13822, 13852, 25337,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35558, 0, 3,
                                                                       34823, 24857, 34844,
                                                                       13852, 13882, 25382,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35621, 0, 3,
                                                                       34844, 24872, 34865,
                                                                       13882, 13912, 25427,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35684, 0, 3,
                                                                       34865, 24887, 34886,
                                                                       13912, 13942, 25472,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35747, 0, 3,
                                                                       34886, 24902, 34907,
                                                                       13942, 13972, 25517,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 35810, 0, 3,
                                                                       34928, 24932, 34991,
                                                                       14032, 14092, 25562,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 35936, 0, 3,
                                                                       34991, 24977, 35054,
                                                                       14092, 14152, 25652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36062, 0, 3,
                                                                       35054, 25022, 35117,
                                                                       14152, 14212, 25742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36188, 0, 3,
                                                                       35117, 25067, 35180,
                                                                       14212, 14272, 25832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36314, 0, 3,
                                                                       35180, 25112, 35243,
                                                                       14272, 14332, 25922,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36440, 0, 3,
                                                                       35243, 25157, 35306,
                                                                       14332, 14392, 26012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36566, 0, 3,
                                                                       35369, 25247, 35432,
                                                                       14512, 14572, 26102,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36692, 0, 3,
                                                                       35432, 25292, 35495,
                                                                       14572, 14632, 26192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36818, 0, 3,
                                                                       35495, 25337, 35558,
                                                                       14632, 14692, 26282,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36944, 0, 3,
                                                                       35558, 25382, 35621,
                                                                       14692, 14752, 26372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37070, 0, 3,
                                                                       35621, 25427, 35684,
                                                                       14752, 14812, 26462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37196, 0, 3,
                                                                       35684, 25472, 35747,
                                                                       14812, 14872, 26552,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37322, 0, 3,
                                                                       35810, 25562, 35936,
                                                                       14992, 15092, 26642,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37532, 0, 3,
                                                                       35936, 25652, 36062,
                                                                       15092, 15192, 26792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37742, 0, 3,
                                                                       36062, 25742, 36188,
                                                                       15192, 15292, 26942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37952, 0, 3,
                                                                       36188, 25832, 36314,
                                                                       15292, 15392, 27092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38162, 0, 3,
                                                                       36314, 25922, 36440,
                                                                       15392, 15492, 27242,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38372, 0, 3,
                                                                       36566, 26102, 36692,
                                                                       15692, 15792, 27392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38582, 0, 3,
                                                                       36692, 26192, 36818,
                                                                       15792, 15892, 27542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38792, 0, 3,
                                                                       36818, 26282, 36944,
                                                                       15892, 15992, 27692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39002, 0, 3,
                                                                       36944, 26372, 37070,
                                                                       15992, 16092, 27842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39212, 0, 3,
                                                                       37070, 26462, 37196,
                                                                       16092, 16192, 27992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 39422, 0, 3,
                                                                       37322, 26642, 37532,
                                                                       16392, 16542, 28142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 39737, 0, 3,
                                                                       37532, 26792, 37742,
                                                                       16542, 16692, 28367,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 40052, 0, 3,
                                                                       37742, 26942, 37952,
                                                                       16692, 16842, 28592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 40367, 0, 3,
                                                                       37952, 27092, 38162,
                                                                       16842, 16992, 28817,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 40682, 0, 3,
                                                                       38372, 27392, 38582,
                                                                       17292, 17442, 29042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 40997, 0, 3,
                                                                       38582, 27542, 38792,
                                                                       17442, 17592, 29267,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 41312, 0, 3,
                                                                       38792, 27692, 39002,
                                                                       17592, 17742, 29492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 41627, 0, 3,
                                                                       39002, 27842, 39212,
                                                                       17742, 17892, 29717,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 41942, 0, 3,
                                                                       39422, 28142, 39737,
                                                                       18192, 18402, 29942,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 42383, 0, 3,
                                                                       39737, 28367, 40052,
                                                                       18402, 18612, 30257,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 42824, 0, 3,
                                                                       40052, 28592, 40367,
                                                                       18612, 18822, 30572,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 43265, 0, 3,
                                                                       40682, 29042, 40997,
                                                                       19242, 19452, 30887,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 43706, 0, 3,
                                                                       40997, 29267, 41312,
                                                                       19452, 19662, 31202,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 44147, 0, 3,
                                                                       41312, 29492, 41627,
                                                                       19662, 19872, 31517,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 44588, 0, 3,
                                                                       41942, 29942, 42383,
                                                                       20292, 20572, 31832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 45176, 0, 3,
                                                                       42383, 30257, 42824,
                                                                       20572, 20852, 32252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 45764, 0, 3,
                                                                       43265, 30887, 43706,
                                                                       21412, 21692, 32672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 46352, 0, 3,
                                                                       43706, 31202, 44147,
                                                                       21692, 21972, 33092,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 46940, 0, 3,
                                                                       44588, 31832, 45176,
                                                                       22532, 22892, 33512,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 47696, 0, 3,
                                                                       45764, 32672, 46352,
                                                                       23612, 23972, 34052,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_g_x(buffer, 48452, 37322, 41942, 1, 21, ncols, beta);

                    simdgeo::geom_g_y(buffer, 48767, 37322, 41942, 1, 21, ncols, beta);

                    simdgeo::geom_g_z(buffer, 49082, 37322, 41942, 1, 21, ncols, beta);

                    simdgeo::geom_g_x(buffer, 49397, 38372, 43265, 1, 21, ncols, beta);

                    simdgeo::geom_g_y(buffer, 49712, 38372, 43265, 1, 21, ncols, beta);

                    simdgeo::geom_g_z(buffer, 50027, 38372, 43265, 1, 21, ncols, beta);

                    simdgeo::geom_h_x(buffer, 50342, 39422, 44588, 1, 21, ncols, beta);

                    simdgeo::geom_h_y(buffer, 50783, 39422, 44588, 1, 21, ncols, beta);

                    simdgeo::geom_h_z(buffer, 51224, 39422, 44588, 1, 21, ncols, beta);

                    simdgeo::geom_h_x(buffer, 51665, 40682, 45764, 1, 21, ncols, beta);

                    simdgeo::geom_h_y(buffer, 52106, 40682, 45764, 1, 21, ncols, beta);

                    simdgeo::geom_h_z(buffer, 52547, 40682, 45764, 1, 21, ncols, beta);

                    simdgeo::geom_i_x(buffer, 52988, 41942, 46940, 1, 21, ncols, beta);

                    simdgeo::geom_i_y(buffer, 53576, 41942, 46940, 1, 21, ncols, beta);

                    simdgeo::geom_i_z(buffer, 54164, 41942, 46940, 1, 21, ncols, beta);

                    simdgeo::geom_i_x(buffer, 54752, 43265, 47696, 1, 21, ncols, beta);

                    simdgeo::geom_i_y(buffer, 55340, 43265, 47696, 1, 21, ncols, beta);

                    simdgeo::geom_i_z(buffer, 55928, 43265, 47696, 1, 21, ncols, beta);

                    simdfunc::contract_primitives(buffer, 56516, 48452, 315, ncols);

                    simdfunc::contract_primitives(buffer, 56996, 48767, 315, ncols);

                    simdfunc::contract_primitives(buffer, 57476, 49082, 315, ncols);

                    simdfunc::contract_primitives(buffer, 57956, 39422, 315, ncols);

                    simdfunc::contract_primitives(buffer, 58436, 49397, 315, ncols);

                    simdfunc::contract_primitives(buffer, 58916, 49712, 315, ncols);

                    simdfunc::contract_primitives(buffer, 59396, 50027, 315, ncols);

                    simdfunc::contract_primitives(buffer, 59876, 40682, 315, ncols);

                    simdfunc::contract_primitives(buffer, 60356, 50342, 441, ncols);

                    simdfunc::contract_primitives(buffer, 61028, 50783, 441, ncols);

                    simdfunc::contract_primitives(buffer, 61700, 51224, 441, ncols);

                    simdfunc::contract_primitives(buffer, 62372, 41942, 441, ncols);

                    simdfunc::contract_primitives(buffer, 63044, 51665, 441, ncols);

                    simdfunc::contract_primitives(buffer, 63716, 52106, 441, ncols);

                    simdfunc::contract_primitives(buffer, 64388, 52547, 441, ncols);

                    simdfunc::contract_primitives(buffer, 65060, 43265, 441, ncols);

                    simdfunc::contract_primitives(buffer, 65732, 52988, 588, ncols);

                    simdfunc::contract_primitives(buffer, 66628, 53576, 588, ncols);

                    simdfunc::contract_primitives(buffer, 67524, 54164, 588, ncols);

                    simdfunc::contract_primitives(buffer, 68420, 54752, 588, ncols);

                    simdfunc::contract_primitives(buffer, 69316, 55340, 588, ncols);

                    simdfunc::contract_primitives(buffer, 70212, 55928, 588, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 56831, 56516, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 57311, 56996, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 57791, 57476, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 58271, 57956, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 58751, 58436, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 59231, 58916, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 59711, 59396, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 60191, 59876, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 60797, 60356, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 61469, 61028, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 62141, 61700, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 62813, 62372, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 63485, 63044, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 64157, 63716, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 64829, 64388, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 65501, 65060, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 66320, 65732, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 67216, 66628, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 68112, 67524, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 69008, 68420, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 69904, 69316, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 70800, 70212, 28, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 71108, 56831, 58271, 60797, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 71603, 57311, 58271, 61469, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 72098, 57791, 58271, 62141, 11,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 72593, 58271, 62813, 11, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 73088, 58751, 60191, 63485, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 73583, 59231, 60191, 64157, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 74078, 59711, 60191, 64829, 11,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 74573, 60191, 65501, 11, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 75068, 60797, 62813, 66320, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 75761, 61469, 62813, 67216, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 76454, 62141, 62813, 68112, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 77147, 63485, 65501, 69008, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 77840, 64157, 65501, 69904, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 78533, 64829, 65501, 70800, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 79226, 71108, 72593, 75068, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 80216, 71603, 72593, 75761, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 81206, 72098, 72593, 76454, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 82196, 73088, 74573, 77147, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 83186, 73583, 74573, 77840, 11,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 84176, 74078, 74573, 78533, 11,
                                          nmax);

        simdtrf::transform_g_inner(buffer, 85166, 82196, 6, 11, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 85166, 99, nmax);

        simdtrf::transform_g_inner(buffer, 85166, 83186, 6, 11, nmax);

        simdtrf::transform_d_outer(values + 495 * nvalues + n * npairs, nvalues, buffer, 85166,
                                   99, nmax);

        simdtrf::transform_g_inner(buffer, 85166, 84176, 6, 11, nmax);

        simdtrf::transform_d_outer(values + 990 * nvalues + n * npairs, nvalues, buffer, 85166,
                                   99, nmax);

        simdtrf::transform_g_inner(buffer, 85166, 79226, 6, 11, nmax);

        simdtrf::transform_d_outer(values + 1485 * nvalues + n * npairs, nvalues, buffer, 85166,
                                   99, nmax);

        simdtrf::transform_g_inner(buffer, 85166, 80216, 6, 11, nmax);

        simdtrf::transform_d_outer(values + 1980 * nvalues + n * npairs, nvalues, buffer, 85166,
                                   99, nmax);

        simdtrf::transform_g_inner(buffer, 85166, 81206, 6, 11, nmax);

        simdtrf::transform_d_outer(values + 2475 * nvalues + n * npairs, nvalues, buffer, 85166,
                                   99, nmax);
    }

    for (size_t m = 0; m < 2970; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
