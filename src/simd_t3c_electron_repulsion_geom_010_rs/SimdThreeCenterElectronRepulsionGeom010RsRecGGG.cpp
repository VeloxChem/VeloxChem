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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGGG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
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
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_ggg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_ggg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 178253, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 4374 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 178253, 72788, 25275, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 13,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 21, 3, 13,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 60, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 63, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 66, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 69, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 72, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 75, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 78, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 81, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 84, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 87, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 90, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 93, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 96, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 99, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 102, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 105, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 108, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 111, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 114, 0, 3, 7, 8,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 120, 0, 3, 8, 9,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 126, 0, 3, 9, 10,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 132, 0, 3, 10, 11,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 138, 0, 3, 11, 12,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 144, 0, 3, 12, 13,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 150, 0, 3, 13, 14,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 156, 0, 3, 14, 15,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 162, 0, 3, 15, 16,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 168, 0, 3, 16, 17,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 174, 0, 3, 17, 18,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 180, 0, 3, 18, 19,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 186, 0, 3, 22, 23,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 192, 0, 3, 23, 24,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 198, 0, 3, 24, 25,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 204, 0, 3, 25, 26,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 210, 0, 3, 26, 27,
                                                                       87, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 216, 0, 3, 27, 28,
                                                                       90, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 222, 0, 3, 28, 29,
                                                                       93, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 228, 0, 3, 29, 30,
                                                                       96, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 234, 0, 3, 30, 31,
                                                                       99, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 240, 0, 3, 31, 32,
                                                                       102, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 246, 0, 3, 32, 33,
                                                                       105, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 252, 0, 3, 33, 34,
                                                                       108, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 36, 39,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 39, 42,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 42, 45,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 45, 48,
                                                                       132, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 48, 51,
                                                                       138, 144, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 51, 54,
                                                                       144, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 54, 57,
                                                                       150, 156, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 57, 60,
                                                                       156, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 60, 63,
                                                                       162, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 63, 66,
                                                                       168, 174, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 66, 69,
                                                                       174, 180, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 75, 78,
                                                                       186, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 78, 81,
                                                                       192, 198, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 81, 84,
                                                                       198, 204, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 84, 87,
                                                                       204, 210, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 87, 90,
                                                                       210, 216, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 90, 93,
                                                                       216, 222, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 93, 96,
                                                                       222, 228, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 96, 99,
                                                                       228, 234, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 99,
                                                                       102, 234, 240, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 102,
                                                                       105, 240, 246, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 105,
                                                                       108, 246, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 478, 0, 3, 114,
                                                                       120, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 493, 0, 3, 120,
                                                                       126, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 508, 0, 3, 126,
                                                                       132, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 523, 0, 3, 132,
                                                                       138, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 538, 0, 3, 138,
                                                                       144, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 553, 0, 3, 144,
                                                                       150, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 568, 0, 3, 150,
                                                                       156, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 583, 0, 3, 156,
                                                                       162, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 598, 0, 3, 162,
                                                                       168, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 613, 0, 3, 168,
                                                                       174, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 628, 0, 3, 186,
                                                                       192, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 643, 0, 3, 192,
                                                                       198, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 658, 0, 3, 198,
                                                                       204, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 673, 0, 3, 204,
                                                                       210, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 688, 0, 3, 210,
                                                                       216, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 703, 0, 3, 216,
                                                                       222, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 718, 0, 3, 222,
                                                                       228, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 733, 0, 3, 228,
                                                                       234, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 748, 0, 3, 234,
                                                                       240, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 763, 0, 3, 240,
                                                                       246, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 778, 0, 3, 258,
                                                                       268, 478, 493, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 799, 0, 3, 268,
                                                                       278, 493, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 820, 0, 3, 278,
                                                                       288, 508, 523, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 841, 0, 3, 288,
                                                                       298, 523, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 862, 0, 3, 298,
                                                                       308, 538, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 883, 0, 3, 308,
                                                                       318, 553, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 904, 0, 3, 318,
                                                                       328, 568, 583, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 925, 0, 3, 328,
                                                                       338, 583, 598, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 946, 0, 3, 338,
                                                                       348, 598, 613, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 967, 0, 3, 368,
                                                                       378, 628, 643, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 988, 0, 3, 378,
                                                                       388, 643, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1009, 0, 3, 388,
                                                                       398, 658, 673, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 398,
                                                                       408, 673, 688, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1051, 0, 3, 408,
                                                                       418, 688, 703, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 418,
                                                                       428, 703, 718, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1093, 0, 3, 428,
                                                                       438, 718, 733, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1114, 0, 3, 438,
                                                                       448, 733, 748, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1135, 0, 3, 448,
                                                                       458, 748, 763, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 478,
                                                                       493, 778, 799, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 493,
                                                                       508, 799, 820, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 508,
                                                                       523, 820, 841, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 523,
                                                                       538, 841, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 538,
                                                                       553, 862, 883, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 553,
                                                                       568, 883, 904, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 568,
                                                                       583, 904, 925, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 583,
                                                                       598, 925, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 628,
                                                                       643, 967, 988, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 643,
                                                                       658, 988, 1009, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 658,
                                                                       673, 1009, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1464, 0, 3, 673,
                                                                       688, 1030, 1051, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 688,
                                                                       703, 1051, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 703,
                                                                       718, 1072, 1093, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 718,
                                                                       733, 1093, 1114, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 733,
                                                                       748, 1114, 1135, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 778,
                                                                       799, 1156, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1640, 0, 3, 799,
                                                                       820, 1184, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1676, 0, 3, 820,
                                                                       841, 1212, 1240, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1712, 0, 3, 841,
                                                                       862, 1240, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1748, 0, 3, 862,
                                                                       883, 1268, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1784, 0, 3, 883,
                                                                       904, 1296, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1820, 0, 3, 904,
                                                                       925, 1324, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 967,
                                                                       988, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1892, 0, 3, 988,
                                                                       1009, 1408, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1928, 0, 3, 1009,
                                                                       1030, 1436, 1464, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1964, 0, 3, 1030,
                                                                       1051, 1464, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2000, 0, 3, 1051,
                                                                       1072, 1492, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2036, 0, 3, 1072,
                                                                       1093, 1520, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2072, 0, 3, 1093,
                                                                       1114, 1548, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 1156,
                                                                       1184, 1604, 1640, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2153, 0, 3, 1184,
                                                                       1212, 1640, 1676, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2198, 0, 3, 1212,
                                                                       1240, 1676, 1712, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2243, 0, 3, 1240,
                                                                       1268, 1712, 1748, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2288, 0, 3, 1268,
                                                                       1296, 1748, 1784, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2333, 0, 3, 1296,
                                                                       1324, 1784, 1820, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2378, 0, 3, 1380,
                                                                       1408, 1856, 1892, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2423, 0, 3, 1408,
                                                                       1436, 1892, 1928, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2468, 0, 3, 1436,
                                                                       1464, 1928, 1964, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2513, 0, 3, 1464,
                                                                       1492, 1964, 2000, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2558, 0, 3, 1492,
                                                                       1520, 2000, 2036, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2603, 0, 3, 1520,
                                                                       1548, 2036, 2072, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 1604,
                                                                       1640, 2108, 2153, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2703, 0, 3, 1640,
                                                                       1676, 2153, 2198, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2758, 0, 3, 1676,
                                                                       1712, 2198, 2243, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2813, 0, 3, 1712,
                                                                       1748, 2243, 2288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2868, 0, 3, 1748,
                                                                       1784, 2288, 2333, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2923, 0, 3, 1856,
                                                                       1892, 2378, 2423, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2978, 0, 3, 1892,
                                                                       1928, 2423, 2468, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3033, 0, 3, 1928,
                                                                       1964, 2468, 2513, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3088, 0, 3, 1964,
                                                                       2000, 2513, 2558, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3143, 0, 3, 2000,
                                                                       2036, 2558, 2603, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3198, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3201, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3204, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3207, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3210, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3213, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3216, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3219, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3222, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3225, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3228, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3231, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3234, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3237, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3240, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3243, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3246, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3249, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3252, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3255, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3258, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3261, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3264, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3267, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3270, 3, 9, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3279, 3, 10, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3288, 3, 11, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3297, 3, 12, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3306, 3, 13, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3315, 3, 14, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3324, 3, 15, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3333, 3, 16, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3342, 3, 17, 66,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3351, 3, 18, 69,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3360, 3, 19, 72,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3369, 3, 24, 81,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3378, 3, 25, 84,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3387, 3, 26, 87,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3396, 3, 27, 90,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3405, 3, 28, 93,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3414, 3, 29, 96,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3423, 3, 30, 99,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3432, 3, 31, 102,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3441, 3, 32, 105,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3450, 3, 33, 108,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3459, 3, 34, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3468, 3, 42, 126,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3486, 3, 45, 132,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3504, 3, 48, 138,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3522, 3, 51, 144,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3540, 3, 54, 150,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3558, 3, 57, 156,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3576, 3, 60, 162,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3594, 3, 63, 168,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3612, 3, 66, 174,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3630, 3, 69, 180,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3648, 3, 81, 198,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3666, 3, 84, 204,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3684, 3, 87, 210,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3702, 3, 90, 216,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3720, 3, 93, 222,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3738, 3, 96, 228,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3756, 3, 99, 234,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3774, 3, 102, 240,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3792, 3, 105, 246,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3810, 3, 108, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3828, 3, 126, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3858, 3, 132, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3888, 3, 138, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3918, 3, 144, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3948, 3, 150, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3978, 3, 156, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4008, 3, 162, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4038, 3, 168, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4068, 3, 174, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4098, 3, 198, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4128, 3, 204, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4158, 3, 210, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4188, 3, 216, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4218, 3, 222, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4248, 3, 228, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4278, 3, 234, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4308, 3, 240, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4338, 3, 246, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4368, 3, 278, 508,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4413, 3, 288, 523,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4458, 3, 298, 538,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4503, 3, 308, 553,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4548, 3, 318, 568,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4593, 3, 328, 583,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4638, 3, 338, 598,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4683, 3, 348, 613,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4728, 3, 388, 658,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4773, 3, 398, 673,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4818, 3, 408, 688,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4863, 3, 418, 703,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4908, 3, 428, 718,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4953, 3, 438, 733,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4998, 3, 448, 748,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5043, 3, 458, 763,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5088, 3, 508, 820,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5151, 3, 523, 841,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5214, 3, 538, 862,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5277, 3, 553, 883,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5340, 3, 568, 904,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5403, 3, 583, 925,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5466, 3, 598, 946,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5529, 3, 658,
                                                                       1009, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5592, 3, 673,
                                                                       1030, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5655, 3, 688,
                                                                       1051, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5718, 3, 703,
                                                                       1072, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5781, 3, 718,
                                                                       1093, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5844, 3, 733,
                                                                       1114, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5907, 3, 748,
                                                                       1135, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5970, 3, 820,
                                                                       1212, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6054, 3, 841,
                                                                       1240, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6138, 3, 862,
                                                                       1268, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6222, 3, 883,
                                                                       1296, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6306, 3, 904,
                                                                       1324, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6390, 3, 925,
                                                                       1352, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6474, 3, 1009,
                                                                       1436, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6558, 3, 1030,
                                                                       1464, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6642, 3, 1051,
                                                                       1492, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6726, 3, 1072,
                                                                       1520, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6810, 3, 1093,
                                                                       1548, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6894, 3, 1114,
                                                                       1576, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6978, 3, 1212,
                                                                       1676, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7086, 3, 1240,
                                                                       1712, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7194, 3, 1268,
                                                                       1748, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7302, 3, 1296,
                                                                       1784, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7410, 3, 1324,
                                                                       1820, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7518, 3, 1436,
                                                                       1928, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7626, 3, 1464,
                                                                       1964, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7734, 3, 1492,
                                                                       2000, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7842, 3, 1520,
                                                                       2036, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7950, 3, 1548,
                                                                       2072, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8058, 3, 1676,
                                                                       2198, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8193, 3, 1712,
                                                                       2243, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8328, 3, 1748,
                                                                       2288, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8463, 3, 1784,
                                                                       2333, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8598, 3, 1928,
                                                                       2468, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8733, 3, 1964,
                                                                       2513, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8868, 3, 2000,
                                                                       2558, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9003, 3, 2036,
                                                                       2603, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9138, 3, 2198,
                                                                       2758, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9303, 3, 2243,
                                                                       2813, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9468, 3, 2288,
                                                                       2868, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9633, 3, 2468,
                                                                       3033, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9798, 3, 2513,
                                                                       3088, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9963, 3, 2558,
                                                                       3143, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10128, 3, 7, 8,
                                                                       3198, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10134, 3, 8, 9,
                                                                       3201, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10140, 3, 9, 10,
                                                                       3204, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10146, 3, 10, 11,
                                                                       3207, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10152, 3, 11, 12,
                                                                       3210, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10158, 3, 12, 13,
                                                                       3213, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10164, 3, 13, 14,
                                                                       3216, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10170, 3, 14, 15,
                                                                       3219, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10176, 3, 15, 16,
                                                                       3222, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10182, 3, 16, 17,
                                                                       3225, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10188, 3, 17, 18,
                                                                       3228, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10194, 3, 18, 19,
                                                                       3231, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10200, 3, 22, 23,
                                                                       3234, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10206, 3, 23, 24,
                                                                       3237, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10212, 3, 24, 25,
                                                                       3240, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10218, 3, 25, 26,
                                                                       3243, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10224, 3, 26, 27,
                                                                       3246, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10230, 3, 27, 28,
                                                                       3249, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10236, 3, 28, 29,
                                                                       3252, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10242, 3, 29, 30,
                                                                       3255, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10248, 3, 30, 31,
                                                                       3258, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10254, 3, 31, 32,
                                                                       3261, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10260, 3, 32, 33,
                                                                       3264, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10266, 3, 33, 34,
                                                                       3267, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10272, 0, 3,
                                                                       10128, 3198, 10134, 3270,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10290, 0, 3,
                                                                       10134, 3201, 10140, 3279,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10308, 0, 3,
                                                                       10140, 3204, 10146, 3288,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10326, 0, 3,
                                                                       10146, 3207, 10152, 3297,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10344, 0, 3,
                                                                       10152, 3210, 10158, 3306,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10362, 0, 3,
                                                                       10158, 3213, 10164, 3315,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10380, 0, 3,
                                                                       10164, 3216, 10170, 3324,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10398, 0, 3,
                                                                       10170, 3219, 10176, 3333,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10416, 0, 3,
                                                                       10176, 3222, 10182, 3342,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10434, 0, 3,
                                                                       10182, 3225, 10188, 3351,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10452, 0, 3,
                                                                       10188, 3228, 10194, 3360,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10470, 0, 3,
                                                                       10200, 3234, 10206, 3369,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10488, 0, 3,
                                                                       10206, 3237, 10212, 3378,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10506, 0, 3,
                                                                       10212, 3240, 10218, 3387,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10524, 0, 3,
                                                                       10218, 3243, 10224, 3396,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10542, 0, 3,
                                                                       10224, 3246, 10230, 3405,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10560, 0, 3,
                                                                       10230, 3249, 10236, 3414,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10578, 0, 3,
                                                                       10236, 3252, 10242, 3423,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10596, 0, 3,
                                                                       10242, 3255, 10248, 3432,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10614, 0, 3,
                                                                       10248, 3258, 10254, 3441,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10632, 0, 3,
                                                                       10254, 3261, 10260, 3450,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10650, 0, 3,
                                                                       10260, 3264, 10266, 3459,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10668, 0, 3,
                                                                       10272, 3270, 10290, 114,
                                                                       120, 3468, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10704, 0, 3,
                                                                       10290, 3279, 10308, 120,
                                                                       126, 3486, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10740, 0, 3,
                                                                       10308, 3288, 10326, 126,
                                                                       132, 3504, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10776, 0, 3,
                                                                       10326, 3297, 10344, 132,
                                                                       138, 3522, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10812, 0, 3,
                                                                       10344, 3306, 10362, 138,
                                                                       144, 3540, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10848, 0, 3,
                                                                       10362, 3315, 10380, 144,
                                                                       150, 3558, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10884, 0, 3,
                                                                       10380, 3324, 10398, 150,
                                                                       156, 3576, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10920, 0, 3,
                                                                       10398, 3333, 10416, 156,
                                                                       162, 3594, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10956, 0, 3,
                                                                       10416, 3342, 10434, 162,
                                                                       168, 3612, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 10992, 0, 3,
                                                                       10434, 3351, 10452, 168,
                                                                       174, 3630, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11028, 0, 3,
                                                                       10470, 3369, 10488, 186,
                                                                       192, 3648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11064, 0, 3,
                                                                       10488, 3378, 10506, 192,
                                                                       198, 3666, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11100, 0, 3,
                                                                       10506, 3387, 10524, 198,
                                                                       204, 3684, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11136, 0, 3,
                                                                       10524, 3396, 10542, 204,
                                                                       210, 3702, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11172, 0, 3,
                                                                       10542, 3405, 10560, 210,
                                                                       216, 3720, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11208, 0, 3,
                                                                       10560, 3414, 10578, 216,
                                                                       222, 3738, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11244, 0, 3,
                                                                       10578, 3423, 10596, 222,
                                                                       228, 3756, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11280, 0, 3,
                                                                       10596, 3432, 10614, 228,
                                                                       234, 3774, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11316, 0, 3,
                                                                       10614, 3441, 10632, 234,
                                                                       240, 3792, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11352, 0, 3,
                                                                       10632, 3450, 10650, 240,
                                                                       246, 3810, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11388, 0, 3,
                                                                       10668, 3468, 10704, 258,
                                                                       268, 3828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11448, 0, 3,
                                                                       10704, 3486, 10740, 268,
                                                                       278, 3858, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11508, 0, 3,
                                                                       10740, 3504, 10776, 278,
                                                                       288, 3888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11568, 0, 3,
                                                                       10776, 3522, 10812, 288,
                                                                       298, 3918, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11628, 0, 3,
                                                                       10812, 3540, 10848, 298,
                                                                       308, 3948, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11688, 0, 3,
                                                                       10848, 3558, 10884, 308,
                                                                       318, 3978, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11748, 0, 3,
                                                                       10884, 3576, 10920, 318,
                                                                       328, 4008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11808, 0, 3,
                                                                       10920, 3594, 10956, 328,
                                                                       338, 4038, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11868, 0, 3,
                                                                       10956, 3612, 10992, 338,
                                                                       348, 4068, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11928, 0, 3,
                                                                       11028, 3648, 11064, 368,
                                                                       378, 4098, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11988, 0, 3,
                                                                       11064, 3666, 11100, 378,
                                                                       388, 4128, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12048, 0, 3,
                                                                       11100, 3684, 11136, 388,
                                                                       398, 4158, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12108, 0, 3,
                                                                       11136, 3702, 11172, 398,
                                                                       408, 4188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12168, 0, 3,
                                                                       11172, 3720, 11208, 408,
                                                                       418, 4218, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12228, 0, 3,
                                                                       11208, 3738, 11244, 418,
                                                                       428, 4248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12288, 0, 3,
                                                                       11244, 3756, 11280, 428,
                                                                       438, 4278, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12348, 0, 3,
                                                                       11280, 3774, 11316, 438,
                                                                       448, 4308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12408, 0, 3,
                                                                       11316, 3792, 11352, 448,
                                                                       458, 4338, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12468, 0, 3,
                                                                       11388, 3828, 11448, 478,
                                                                       493, 4368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12558, 0, 3,
                                                                       11448, 3858, 11508, 493,
                                                                       508, 4413, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12648, 0, 3,
                                                                       11508, 3888, 11568, 508,
                                                                       523, 4458, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12738, 0, 3,
                                                                       11568, 3918, 11628, 523,
                                                                       538, 4503, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12828, 0, 3,
                                                                       11628, 3948, 11688, 538,
                                                                       553, 4548, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12918, 0, 3,
                                                                       11688, 3978, 11748, 553,
                                                                       568, 4593, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13008, 0, 3,
                                                                       11748, 4008, 11808, 568,
                                                                       583, 4638, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13098, 0, 3,
                                                                       11808, 4038, 11868, 583,
                                                                       598, 4683, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13188, 0, 3,
                                                                       11928, 4098, 11988, 628,
                                                                       643, 4728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13278, 0, 3,
                                                                       11988, 4128, 12048, 643,
                                                                       658, 4773, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13368, 0, 3,
                                                                       12048, 4158, 12108, 658,
                                                                       673, 4818, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13458, 0, 3,
                                                                       12108, 4188, 12168, 673,
                                                                       688, 4863, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13548, 0, 3,
                                                                       12168, 4218, 12228, 688,
                                                                       703, 4908, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13638, 0, 3,
                                                                       12228, 4248, 12288, 703,
                                                                       718, 4953, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13728, 0, 3,
                                                                       12288, 4278, 12348, 718,
                                                                       733, 4998, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13818, 0, 3,
                                                                       12348, 4308, 12408, 733,
                                                                       748, 5043, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13908, 0, 3,
                                                                       12468, 4368, 12558, 778,
                                                                       799, 5088, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14034, 0, 3,
                                                                       12558, 4413, 12648, 799,
                                                                       820, 5151, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14160, 0, 3,
                                                                       12648, 4458, 12738, 820,
                                                                       841, 5214, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14286, 0, 3,
                                                                       12738, 4503, 12828, 841,
                                                                       862, 5277, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14412, 0, 3,
                                                                       12828, 4548, 12918, 862,
                                                                       883, 5340, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14538, 0, 3,
                                                                       12918, 4593, 13008, 883,
                                                                       904, 5403, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14664, 0, 3,
                                                                       13008, 4638, 13098, 904,
                                                                       925, 5466, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14790, 0, 3,
                                                                       13188, 4728, 13278, 967,
                                                                       988, 5529, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14916, 0, 3,
                                                                       13278, 4773, 13368, 988,
                                                                       1009, 5592, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15042, 0, 3,
                                                                       13368, 4818, 13458, 1009,
                                                                       1030, 5655, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15168, 0, 3,
                                                                       13458, 4863, 13548, 1030,
                                                                       1051, 5718, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15294, 0, 3,
                                                                       13548, 4908, 13638, 1051,
                                                                       1072, 5781, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15420, 0, 3,
                                                                       13638, 4953, 13728, 1072,
                                                                       1093, 5844, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15546, 0, 3,
                                                                       13728, 4998, 13818, 1093,
                                                                       1114, 5907, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15672, 0, 3,
                                                                       13908, 5088, 14034, 1156,
                                                                       1184, 5970, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15840, 0, 3,
                                                                       14034, 5151, 14160, 1184,
                                                                       1212, 6054, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16008, 0, 3,
                                                                       14160, 5214, 14286, 1212,
                                                                       1240, 6138, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16176, 0, 3,
                                                                       14286, 5277, 14412, 1240,
                                                                       1268, 6222, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16344, 0, 3,
                                                                       14412, 5340, 14538, 1268,
                                                                       1296, 6306, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16512, 0, 3,
                                                                       14538, 5403, 14664, 1296,
                                                                       1324, 6390, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16680, 0, 3,
                                                                       14790, 5529, 14916, 1380,
                                                                       1408, 6474, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16848, 0, 3,
                                                                       14916, 5592, 15042, 1408,
                                                                       1436, 6558, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17016, 0, 3,
                                                                       15042, 5655, 15168, 1436,
                                                                       1464, 6642, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17184, 0, 3,
                                                                       15168, 5718, 15294, 1464,
                                                                       1492, 6726, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17352, 0, 3,
                                                                       15294, 5781, 15420, 1492,
                                                                       1520, 6810, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17520, 0, 3,
                                                                       15420, 5844, 15546, 1520,
                                                                       1548, 6894, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17688, 0, 3,
                                                                       15672, 5970, 15840, 1604,
                                                                       1640, 6978, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17904, 0, 3,
                                                                       15840, 6054, 16008, 1640,
                                                                       1676, 7086, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 18120, 0, 3,
                                                                       16008, 6138, 16176, 1676,
                                                                       1712, 7194, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 18336, 0, 3,
                                                                       16176, 6222, 16344, 1712,
                                                                       1748, 7302, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 18552, 0, 3,
                                                                       16344, 6306, 16512, 1748,
                                                                       1784, 7410, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 18768, 0, 3,
                                                                       16680, 6474, 16848, 1856,
                                                                       1892, 7518, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 18984, 0, 3,
                                                                       16848, 6558, 17016, 1892,
                                                                       1928, 7626, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19200, 0, 3,
                                                                       17016, 6642, 17184, 1928,
                                                                       1964, 7734, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19416, 0, 3,
                                                                       17184, 6726, 17352, 1964,
                                                                       2000, 7842, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19632, 0, 3,
                                                                       17352, 6810, 17520, 2000,
                                                                       2036, 7950, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 19848, 0, 3,
                                                                       17688, 6978, 17904, 2108,
                                                                       2153, 8058, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 20118, 0, 3,
                                                                       17904, 7086, 18120, 2153,
                                                                       2198, 8193, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 20388, 0, 3,
                                                                       18120, 7194, 18336, 2198,
                                                                       2243, 8328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 20658, 0, 3,
                                                                       18336, 7302, 18552, 2243,
                                                                       2288, 8463, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 20928, 0, 3,
                                                                       18768, 7518, 18984, 2378,
                                                                       2423, 8598, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 21198, 0, 3,
                                                                       18984, 7626, 19200, 2423,
                                                                       2468, 8733, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 21468, 0, 3,
                                                                       19200, 7734, 19416, 2468,
                                                                       2513, 8868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 21738, 0, 3,
                                                                       19416, 7842, 19632, 2513,
                                                                       2558, 9003, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 22008, 0, 3,
                                                                       19848, 8058, 20118, 2648,
                                                                       2703, 9138, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 22338, 0, 3,
                                                                       20118, 8193, 20388, 2703,
                                                                       2758, 9303, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 22668, 0, 3,
                                                                       20388, 8328, 20658, 2758,
                                                                       2813, 9468, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 22998, 0, 3,
                                                                       20928, 8598, 21198, 2923,
                                                                       2978, 9633, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 23328, 0, 3,
                                                                       21198, 8733, 21468, 2978,
                                                                       3033, 9798, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 23658, 0, 3,
                                                                       21468, 8868, 21738, 3033,
                                                                       3088, 9963, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23988, 3, 3198,
                                                                       3201, 10140, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23998, 3, 3201,
                                                                       3204, 10146, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24008, 3, 3204,
                                                                       3207, 10152, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24018, 3, 3207,
                                                                       3210, 10158, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24028, 3, 3210,
                                                                       3213, 10164, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24038, 3, 3213,
                                                                       3216, 10170, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24048, 3, 3216,
                                                                       3219, 10176, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24058, 3, 3219,
                                                                       3222, 10182, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24068, 3, 3222,
                                                                       3225, 10188, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24078, 3, 3225,
                                                                       3228, 10194, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24088, 3, 3234,
                                                                       3237, 10212, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24098, 3, 3237,
                                                                       3240, 10218, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24108, 3, 3240,
                                                                       3243, 10224, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24118, 3, 3243,
                                                                       3246, 10230, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24128, 3, 3246,
                                                                       3249, 10236, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24138, 3, 3249,
                                                                       3252, 10242, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24148, 3, 3252,
                                                                       3255, 10248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24158, 3, 3255,
                                                                       3258, 10254, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24168, 3, 3258,
                                                                       3261, 10260, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24178, 3, 3261,
                                                                       3264, 10266, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24188, 0, 3,
                                                                       23988, 10140, 23998,
                                                                       10308, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24218, 0, 3,
                                                                       23998, 10146, 24008,
                                                                       10326, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24248, 0, 3,
                                                                       24008, 10152, 24018,
                                                                       10344, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24278, 0, 3,
                                                                       24018, 10158, 24028,
                                                                       10362, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24308, 0, 3,
                                                                       24028, 10164, 24038,
                                                                       10380, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24338, 0, 3,
                                                                       24038, 10170, 24048,
                                                                       10398, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24368, 0, 3,
                                                                       24048, 10176, 24058,
                                                                       10416, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24398, 0, 3,
                                                                       24058, 10182, 24068,
                                                                       10434, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24428, 0, 3,
                                                                       24068, 10188, 24078,
                                                                       10452, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24458, 0, 3,
                                                                       24088, 10212, 24098,
                                                                       10506, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24488, 0, 3,
                                                                       24098, 10218, 24108,
                                                                       10524, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24518, 0, 3,
                                                                       24108, 10224, 24118,
                                                                       10542, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24548, 0, 3,
                                                                       24118, 10230, 24128,
                                                                       10560, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24578, 0, 3,
                                                                       24128, 10236, 24138,
                                                                       10578, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24608, 0, 3,
                                                                       24138, 10242, 24148,
                                                                       10596, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24638, 0, 3,
                                                                       24148, 10248, 24158,
                                                                       10614, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24668, 0, 3,
                                                                       24158, 10254, 24168,
                                                                       10632, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24698, 0, 3,
                                                                       24168, 10260, 24178,
                                                                       10650, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24728, 0, 3,
                                                                       24188, 10308, 24218, 3468,
                                                                       3486, 10740, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24788, 0, 3,
                                                                       24218, 10326, 24248, 3486,
                                                                       3504, 10776, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24848, 0, 3,
                                                                       24248, 10344, 24278, 3504,
                                                                       3522, 10812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24908, 0, 3,
                                                                       24278, 10362, 24308, 3522,
                                                                       3540, 10848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24968, 0, 3,
                                                                       24308, 10380, 24338, 3540,
                                                                       3558, 10884, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25028, 0, 3,
                                                                       24338, 10398, 24368, 3558,
                                                                       3576, 10920, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25088, 0, 3,
                                                                       24368, 10416, 24398, 3576,
                                                                       3594, 10956, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25148, 0, 3,
                                                                       24398, 10434, 24428, 3594,
                                                                       3612, 10992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25208, 0, 3,
                                                                       24458, 10506, 24488, 3648,
                                                                       3666, 11100, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25268, 0, 3,
                                                                       24488, 10524, 24518, 3666,
                                                                       3684, 11136, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25328, 0, 3,
                                                                       24518, 10542, 24548, 3684,
                                                                       3702, 11172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25388, 0, 3,
                                                                       24548, 10560, 24578, 3702,
                                                                       3720, 11208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25448, 0, 3,
                                                                       24578, 10578, 24608, 3720,
                                                                       3738, 11244, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25508, 0, 3,
                                                                       24608, 10596, 24638, 3738,
                                                                       3756, 11280, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25568, 0, 3,
                                                                       24638, 10614, 24668, 3756,
                                                                       3774, 11316, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25628, 0, 3,
                                                                       24668, 10632, 24698, 3774,
                                                                       3792, 11352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25688, 0, 3,
                                                                       24728, 10740, 24788, 3828,
                                                                       3858, 11508, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25788, 0, 3,
                                                                       24788, 10776, 24848, 3858,
                                                                       3888, 11568, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25888, 0, 3,
                                                                       24848, 10812, 24908, 3888,
                                                                       3918, 11628, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25988, 0, 3,
                                                                       24908, 10848, 24968, 3918,
                                                                       3948, 11688, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26088, 0, 3,
                                                                       24968, 10884, 25028, 3948,
                                                                       3978, 11748, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26188, 0, 3,
                                                                       25028, 10920, 25088, 3978,
                                                                       4008, 11808, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26288, 0, 3,
                                                                       25088, 10956, 25148, 4008,
                                                                       4038, 11868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26388, 0, 3,
                                                                       25208, 11100, 25268, 4098,
                                                                       4128, 12048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26488, 0, 3,
                                                                       25268, 11136, 25328, 4128,
                                                                       4158, 12108, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26588, 0, 3,
                                                                       25328, 11172, 25388, 4158,
                                                                       4188, 12168, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26688, 0, 3,
                                                                       25388, 11208, 25448, 4188,
                                                                       4218, 12228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26788, 0, 3,
                                                                       25448, 11244, 25508, 4218,
                                                                       4248, 12288, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26888, 0, 3,
                                                                       25508, 11280, 25568, 4248,
                                                                       4278, 12348, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26988, 0, 3,
                                                                       25568, 11316, 25628, 4278,
                                                                       4308, 12408, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27088, 0, 3,
                                                                       25688, 11508, 25788, 4368,
                                                                       4413, 12648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27238, 0, 3,
                                                                       25788, 11568, 25888, 4413,
                                                                       4458, 12738, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27388, 0, 3,
                                                                       25888, 11628, 25988, 4458,
                                                                       4503, 12828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27538, 0, 3,
                                                                       25988, 11688, 26088, 4503,
                                                                       4548, 12918, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27688, 0, 3,
                                                                       26088, 11748, 26188, 4548,
                                                                       4593, 13008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27838, 0, 3,
                                                                       26188, 11808, 26288, 4593,
                                                                       4638, 13098, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27988, 0, 3,
                                                                       26388, 12048, 26488, 4728,
                                                                       4773, 13368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 28138, 0, 3,
                                                                       26488, 12108, 26588, 4773,
                                                                       4818, 13458, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 28288, 0, 3,
                                                                       26588, 12168, 26688, 4818,
                                                                       4863, 13548, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 28438, 0, 3,
                                                                       26688, 12228, 26788, 4863,
                                                                       4908, 13638, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 28588, 0, 3,
                                                                       26788, 12288, 26888, 4908,
                                                                       4953, 13728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 28738, 0, 3,
                                                                       26888, 12348, 26988, 4953,
                                                                       4998, 13818, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28888, 0, 3,
                                                                       27088, 12648, 27238, 5088,
                                                                       5151, 14160, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29098, 0, 3,
                                                                       27238, 12738, 27388, 5151,
                                                                       5214, 14286, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29308, 0, 3,
                                                                       27388, 12828, 27538, 5214,
                                                                       5277, 14412, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29518, 0, 3,
                                                                       27538, 12918, 27688, 5277,
                                                                       5340, 14538, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29728, 0, 3,
                                                                       27688, 13008, 27838, 5340,
                                                                       5403, 14664, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29938, 0, 3,
                                                                       27988, 13368, 28138, 5529,
                                                                       5592, 15042, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 30148, 0, 3,
                                                                       28138, 13458, 28288, 5592,
                                                                       5655, 15168, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 30358, 0, 3,
                                                                       28288, 13548, 28438, 5655,
                                                                       5718, 15294, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 30568, 0, 3,
                                                                       28438, 13638, 28588, 5718,
                                                                       5781, 15420, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 30778, 0, 3,
                                                                       28588, 13728, 28738, 5781,
                                                                       5844, 15546, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30988, 0, 3,
                                                                       28888, 14160, 29098, 5970,
                                                                       6054, 16008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31268, 0, 3,
                                                                       29098, 14286, 29308, 6054,
                                                                       6138, 16176, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31548, 0, 3,
                                                                       29308, 14412, 29518, 6138,
                                                                       6222, 16344, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31828, 0, 3,
                                                                       29518, 14538, 29728, 6222,
                                                                       6306, 16512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 32108, 0, 3,
                                                                       29938, 15042, 30148, 6474,
                                                                       6558, 17016, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 32388, 0, 3,
                                                                       30148, 15168, 30358, 6558,
                                                                       6642, 17184, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 32668, 0, 3,
                                                                       30358, 15294, 30568, 6642,
                                                                       6726, 17352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 32948, 0, 3,
                                                                       30568, 15420, 30778, 6726,
                                                                       6810, 17520, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33228, 0, 3,
                                                                       30988, 16008, 31268, 6978,
                                                                       7086, 18120, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33588, 0, 3,
                                                                       31268, 16176, 31548, 7086,
                                                                       7194, 18336, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33948, 0, 3,
                                                                       31548, 16344, 31828, 7194,
                                                                       7302, 18552, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 34308, 0, 3,
                                                                       32108, 17016, 32388, 7518,
                                                                       7626, 19200, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 34668, 0, 3,
                                                                       32388, 17184, 32668, 7626,
                                                                       7734, 19416, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 35028, 0, 3,
                                                                       32668, 17352, 32948, 7734,
                                                                       7842, 19632, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 35388, 0, 3,
                                                                       33228, 18120, 33588, 8058,
                                                                       8193, 20388, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 35838, 0, 3,
                                                                       33588, 18336, 33948, 8193,
                                                                       8328, 20658, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 36288, 0, 3,
                                                                       34308, 19200, 34668, 8598,
                                                                       8733, 21468, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 36738, 0, 3,
                                                                       34668, 19416, 35028, 8733,
                                                                       8868, 21738, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 37188, 0, 3,
                                                                       35388, 20388, 35838, 9138,
                                                                       9303, 22668, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 37738, 0, 3,
                                                                       36288, 21468, 36738, 9633,
                                                                       9798, 23658, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38288, 3, 10128,
                                                                       10134, 23988, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38303, 3, 10134,
                                                                       10140, 23998, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38318, 3, 10140,
                                                                       10146, 24008, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38333, 3, 10146,
                                                                       10152, 24018, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38348, 3, 10152,
                                                                       10158, 24028, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38363, 3, 10158,
                                                                       10164, 24038, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38378, 3, 10164,
                                                                       10170, 24048, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38393, 3, 10170,
                                                                       10176, 24058, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38408, 3, 10176,
                                                                       10182, 24068, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38423, 3, 10182,
                                                                       10188, 24078, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38438, 3, 10200,
                                                                       10206, 24088, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38453, 3, 10206,
                                                                       10212, 24098, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38468, 3, 10212,
                                                                       10218, 24108, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38483, 3, 10218,
                                                                       10224, 24118, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38498, 3, 10224,
                                                                       10230, 24128, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38513, 3, 10230,
                                                                       10236, 24138, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38528, 3, 10236,
                                                                       10242, 24148, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38543, 3, 10242,
                                                                       10248, 24158, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38558, 3, 10248,
                                                                       10254, 24168, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 38573, 3, 10254,
                                                                       10260, 24178, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38588, 0, 3,
                                                                       38288, 23988, 38303,
                                                                       10272, 10290, 24188,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38633, 0, 3,
                                                                       38303, 23998, 38318,
                                                                       10290, 10308, 24218,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38678, 0, 3,
                                                                       38318, 24008, 38333,
                                                                       10308, 10326, 24248,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38723, 0, 3,
                                                                       38333, 24018, 38348,
                                                                       10326, 10344, 24278,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38768, 0, 3,
                                                                       38348, 24028, 38363,
                                                                       10344, 10362, 24308,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38813, 0, 3,
                                                                       38363, 24038, 38378,
                                                                       10362, 10380, 24338,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38858, 0, 3,
                                                                       38378, 24048, 38393,
                                                                       10380, 10398, 24368,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38903, 0, 3,
                                                                       38393, 24058, 38408,
                                                                       10398, 10416, 24398,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38948, 0, 3,
                                                                       38408, 24068, 38423,
                                                                       10416, 10434, 24428,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38993, 0, 3,
                                                                       38438, 24088, 38453,
                                                                       10470, 10488, 24458,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 39038, 0, 3,
                                                                       38453, 24098, 38468,
                                                                       10488, 10506, 24488,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 39083, 0, 3,
                                                                       38468, 24108, 38483,
                                                                       10506, 10524, 24518,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 39128, 0, 3,
                                                                       38483, 24118, 38498,
                                                                       10524, 10542, 24548,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 39173, 0, 3,
                                                                       38498, 24128, 38513,
                                                                       10542, 10560, 24578,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 39218, 0, 3,
                                                                       38513, 24138, 38528,
                                                                       10560, 10578, 24608,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 39263, 0, 3,
                                                                       38528, 24148, 38543,
                                                                       10578, 10596, 24638,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 39308, 0, 3,
                                                                       38543, 24158, 38558,
                                                                       10596, 10614, 24668,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 39353, 0, 3,
                                                                       38558, 24168, 38573,
                                                                       10614, 10632, 24698,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39398, 0, 3,
                                                                       38588, 24188, 38633,
                                                                       10668, 10704, 24728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39488, 0, 3,
                                                                       38633, 24218, 38678,
                                                                       10704, 10740, 24788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39578, 0, 3,
                                                                       38678, 24248, 38723,
                                                                       10740, 10776, 24848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39668, 0, 3,
                                                                       38723, 24278, 38768,
                                                                       10776, 10812, 24908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39758, 0, 3,
                                                                       38768, 24308, 38813,
                                                                       10812, 10848, 24968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39848, 0, 3,
                                                                       38813, 24338, 38858,
                                                                       10848, 10884, 25028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39938, 0, 3,
                                                                       38858, 24368, 38903,
                                                                       10884, 10920, 25088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40028, 0, 3,
                                                                       38903, 24398, 38948,
                                                                       10920, 10956, 25148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40118, 0, 3,
                                                                       38993, 24458, 39038,
                                                                       11028, 11064, 25208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40208, 0, 3,
                                                                       39038, 24488, 39083,
                                                                       11064, 11100, 25268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40298, 0, 3,
                                                                       39083, 24518, 39128,
                                                                       11100, 11136, 25328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40388, 0, 3,
                                                                       39128, 24548, 39173,
                                                                       11136, 11172, 25388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40478, 0, 3,
                                                                       39173, 24578, 39218,
                                                                       11172, 11208, 25448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40568, 0, 3,
                                                                       39218, 24608, 39263,
                                                                       11208, 11244, 25508,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40658, 0, 3,
                                                                       39263, 24638, 39308,
                                                                       11244, 11280, 25568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 40748, 0, 3,
                                                                       39308, 24668, 39353,
                                                                       11280, 11316, 25628,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40838, 0, 3,
                                                                       39398, 24728, 39488,
                                                                       11388, 11448, 25688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40988, 0, 3,
                                                                       39488, 24788, 39578,
                                                                       11448, 11508, 25788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 41138, 0, 3,
                                                                       39578, 24848, 39668,
                                                                       11508, 11568, 25888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 41288, 0, 3,
                                                                       39668, 24908, 39758,
                                                                       11568, 11628, 25988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 41438, 0, 3,
                                                                       39758, 24968, 39848,
                                                                       11628, 11688, 26088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 41588, 0, 3,
                                                                       39848, 25028, 39938,
                                                                       11688, 11748, 26188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 41738, 0, 3,
                                                                       39938, 25088, 40028,
                                                                       11748, 11808, 26288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 41888, 0, 3,
                                                                       40118, 25208, 40208,
                                                                       11928, 11988, 26388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 42038, 0, 3,
                                                                       40208, 25268, 40298,
                                                                       11988, 12048, 26488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 42188, 0, 3,
                                                                       40298, 25328, 40388,
                                                                       12048, 12108, 26588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 42338, 0, 3,
                                                                       40388, 25388, 40478,
                                                                       12108, 12168, 26688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 42488, 0, 3,
                                                                       40478, 25448, 40568,
                                                                       12168, 12228, 26788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 42638, 0, 3,
                                                                       40568, 25508, 40658,
                                                                       12228, 12288, 26888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 42788, 0, 3,
                                                                       40658, 25568, 40748,
                                                                       12288, 12348, 26988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 42938, 0, 3,
                                                                       40838, 25688, 40988,
                                                                       12468, 12558, 27088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 43163, 0, 3,
                                                                       40988, 25788, 41138,
                                                                       12558, 12648, 27238,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 43388, 0, 3,
                                                                       41138, 25888, 41288,
                                                                       12648, 12738, 27388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 43613, 0, 3,
                                                                       41288, 25988, 41438,
                                                                       12738, 12828, 27538,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 43838, 0, 3,
                                                                       41438, 26088, 41588,
                                                                       12828, 12918, 27688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 44063, 0, 3,
                                                                       41588, 26188, 41738,
                                                                       12918, 13008, 27838,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 44288, 0, 3,
                                                                       41888, 26388, 42038,
                                                                       13188, 13278, 27988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 44513, 0, 3,
                                                                       42038, 26488, 42188,
                                                                       13278, 13368, 28138,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 44738, 0, 3,
                                                                       42188, 26588, 42338,
                                                                       13368, 13458, 28288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 44963, 0, 3,
                                                                       42338, 26688, 42488,
                                                                       13458, 13548, 28438,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 45188, 0, 3,
                                                                       42488, 26788, 42638,
                                                                       13548, 13638, 28588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 45413, 0, 3,
                                                                       42638, 26888, 42788,
                                                                       13638, 13728, 28738,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 45638, 0, 3,
                                                                       42938, 27088, 43163,
                                                                       13908, 14034, 28888,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 45953, 0, 3,
                                                                       43163, 27238, 43388,
                                                                       14034, 14160, 29098,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 46268, 0, 3,
                                                                       43388, 27388, 43613,
                                                                       14160, 14286, 29308,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 46583, 0, 3,
                                                                       43613, 27538, 43838,
                                                                       14286, 14412, 29518,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 46898, 0, 3,
                                                                       43838, 27688, 44063,
                                                                       14412, 14538, 29728,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 47213, 0, 3,
                                                                       44288, 27988, 44513,
                                                                       14790, 14916, 29938,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 47528, 0, 3,
                                                                       44513, 28138, 44738,
                                                                       14916, 15042, 30148,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 47843, 0, 3,
                                                                       44738, 28288, 44963,
                                                                       15042, 15168, 30358,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 48158, 0, 3,
                                                                       44963, 28438, 45188,
                                                                       15168, 15294, 30568,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 48473, 0, 3,
                                                                       45188, 28588, 45413,
                                                                       15294, 15420, 30778,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 48788, 0, 3,
                                                                       45638, 28888, 45953,
                                                                       15672, 15840, 30988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 49208, 0, 3,
                                                                       45953, 29098, 46268,
                                                                       15840, 16008, 31268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 49628, 0, 3,
                                                                       46268, 29308, 46583,
                                                                       16008, 16176, 31548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 50048, 0, 3,
                                                                       46583, 29518, 46898,
                                                                       16176, 16344, 31828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 50468, 0, 3,
                                                                       47213, 29938, 47528,
                                                                       16680, 16848, 32108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 50888, 0, 3,
                                                                       47528, 30148, 47843,
                                                                       16848, 17016, 32388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 51308, 0, 3,
                                                                       47843, 30358, 48158,
                                                                       17016, 17184, 32668,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 51728, 0, 3,
                                                                       48158, 30568, 48473,
                                                                       17184, 17352, 32948,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 52148, 0, 3,
                                                                       48788, 30988, 49208,
                                                                       17688, 17904, 33228,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 52688, 0, 3,
                                                                       49208, 31268, 49628,
                                                                       17904, 18120, 33588,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 53228, 0, 3,
                                                                       49628, 31548, 50048,
                                                                       18120, 18336, 33948,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 53768, 0, 3,
                                                                       50468, 32108, 50888,
                                                                       18768, 18984, 34308,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 54308, 0, 3,
                                                                       50888, 32388, 51308,
                                                                       18984, 19200, 34668,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 54848, 0, 3,
                                                                       51308, 32668, 51728,
                                                                       19200, 19416, 35028,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 55388, 0, 3,
                                                                       52148, 33228, 52688,
                                                                       19848, 20118, 35388,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 56063, 0, 3,
                                                                       52688, 33588, 53228,
                                                                       20118, 20388, 35838,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 56738, 0, 3,
                                                                       53768, 34308, 54308,
                                                                       20928, 21198, 36288,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 57413, 0, 3,
                                                                       54308, 34668, 54848,
                                                                       21198, 21468, 36738,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 58088, 0, 3,
                                                                       55388, 35388, 56063,
                                                                       22008, 22338, 37188,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 58913, 0, 3,
                                                                       56738, 36288, 57413,
                                                                       22998, 23328, 37738,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_g_x(buffer, 59738, 40838, 45638, 1, 15, ncols, beta);

                    simdgeo::geom_g_y(buffer, 59963, 40838, 45638, 1, 15, ncols, beta);

                    simdgeo::geom_g_z(buffer, 60188, 40838, 45638, 1, 15, ncols, beta);

                    simdgeo::geom_g_x(buffer, 60413, 41888, 47213, 1, 15, ncols, beta);

                    simdgeo::geom_g_y(buffer, 60638, 41888, 47213, 1, 15, ncols, beta);

                    simdgeo::geom_g_z(buffer, 60863, 41888, 47213, 1, 15, ncols, beta);

                    simdgeo::geom_h_x(buffer, 61088, 42938, 48788, 1, 15, ncols, beta);

                    simdgeo::geom_h_y(buffer, 61403, 42938, 48788, 1, 15, ncols, beta);

                    simdgeo::geom_h_z(buffer, 61718, 42938, 48788, 1, 15, ncols, beta);

                    simdgeo::geom_h_x(buffer, 62033, 44288, 50468, 1, 15, ncols, beta);

                    simdgeo::geom_h_y(buffer, 62348, 44288, 50468, 1, 15, ncols, beta);

                    simdgeo::geom_h_z(buffer, 62663, 44288, 50468, 1, 15, ncols, beta);

                    simdgeo::geom_i_x(buffer, 62978, 45638, 52148, 1, 15, ncols, beta);

                    simdgeo::geom_i_y(buffer, 63398, 45638, 52148, 1, 15, ncols, beta);

                    simdgeo::geom_i_z(buffer, 63818, 45638, 52148, 1, 15, ncols, beta);

                    simdgeo::geom_i_x(buffer, 64238, 47213, 53768, 1, 15, ncols, beta);

                    simdgeo::geom_i_y(buffer, 64658, 47213, 53768, 1, 15, ncols, beta);

                    simdgeo::geom_i_z(buffer, 65078, 47213, 53768, 1, 15, ncols, beta);

                    simdgeo::geom_k_x(buffer, 65498, 48788, 55388, 1, 15, ncols, beta);

                    simdgeo::geom_k_y(buffer, 66038, 48788, 55388, 1, 15, ncols, beta);

                    simdgeo::geom_k_z(buffer, 66578, 48788, 55388, 1, 15, ncols, beta);

                    simdgeo::geom_k_x(buffer, 67118, 50468, 56738, 1, 15, ncols, beta);

                    simdgeo::geom_k_y(buffer, 67658, 50468, 56738, 1, 15, ncols, beta);

                    simdgeo::geom_k_z(buffer, 68198, 50468, 56738, 1, 15, ncols, beta);

                    simdgeo::geom_l_x(buffer, 68738, 52148, 58088, 1, 15, ncols, beta);

                    simdgeo::geom_l_y(buffer, 69413, 52148, 58088, 1, 15, ncols, beta);

                    simdgeo::geom_l_z(buffer, 70088, 52148, 58088, 1, 15, ncols, beta);

                    simdgeo::geom_l_x(buffer, 70763, 53768, 58913, 1, 15, ncols, beta);

                    simdgeo::geom_l_y(buffer, 71438, 53768, 58913, 1, 15, ncols, beta);

                    simdgeo::geom_l_z(buffer, 72113, 53768, 58913, 1, 15, ncols, beta);

                    simdfunc::contract_primitives(buffer, 72788, 59738, 225, ncols);

                    simdfunc::contract_primitives(buffer, 73148, 59963, 225, ncols);

                    simdfunc::contract_primitives(buffer, 73508, 60188, 225, ncols);

                    simdfunc::contract_primitives(buffer, 73868, 42938, 225, ncols);

                    simdfunc::contract_primitives(buffer, 74228, 60413, 225, ncols);

                    simdfunc::contract_primitives(buffer, 74588, 60638, 225, ncols);

                    simdfunc::contract_primitives(buffer, 74948, 60863, 225, ncols);

                    simdfunc::contract_primitives(buffer, 75308, 44288, 225, ncols);

                    simdfunc::contract_primitives(buffer, 75668, 61088, 315, ncols);

                    simdfunc::contract_primitives(buffer, 76172, 61403, 315, ncols);

                    simdfunc::contract_primitives(buffer, 76676, 61718, 315, ncols);

                    simdfunc::contract_primitives(buffer, 77180, 45638, 315, ncols);

                    simdfunc::contract_primitives(buffer, 77684, 62033, 315, ncols);

                    simdfunc::contract_primitives(buffer, 78188, 62348, 315, ncols);

                    simdfunc::contract_primitives(buffer, 78692, 62663, 315, ncols);

                    simdfunc::contract_primitives(buffer, 79196, 47213, 315, ncols);

                    simdfunc::contract_primitives(buffer, 79700, 62978, 420, ncols);

                    simdfunc::contract_primitives(buffer, 80372, 63398, 420, ncols);

                    simdfunc::contract_primitives(buffer, 81044, 63818, 420, ncols);

                    simdfunc::contract_primitives(buffer, 81716, 48788, 420, ncols);

                    simdfunc::contract_primitives(buffer, 82388, 64238, 420, ncols);

                    simdfunc::contract_primitives(buffer, 83060, 64658, 420, ncols);

                    simdfunc::contract_primitives(buffer, 83732, 65078, 420, ncols);

                    simdfunc::contract_primitives(buffer, 84404, 50468, 420, ncols);

                    simdfunc::contract_primitives(buffer, 85076, 65498, 540, ncols);

                    simdfunc::contract_primitives(buffer, 85940, 66038, 540, ncols);

                    simdfunc::contract_primitives(buffer, 86804, 66578, 540, ncols);

                    simdfunc::contract_primitives(buffer, 87668, 52148, 540, ncols);

                    simdfunc::contract_primitives(buffer, 88532, 67118, 540, ncols);

                    simdfunc::contract_primitives(buffer, 89396, 67658, 540, ncols);

                    simdfunc::contract_primitives(buffer, 90260, 68198, 540, ncols);

                    simdfunc::contract_primitives(buffer, 91124, 53768, 540, ncols);

                    simdfunc::contract_primitives(buffer, 91988, 68738, 675, ncols);

                    simdfunc::contract_primitives(buffer, 93068, 69413, 675, ncols);

                    simdfunc::contract_primitives(buffer, 94148, 70088, 675, ncols);

                    simdfunc::contract_primitives(buffer, 95228, 70763, 675, ncols);

                    simdfunc::contract_primitives(buffer, 96308, 71438, 675, ncols);

                    simdfunc::contract_primitives(buffer, 97388, 72113, 675, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 73013, 72788, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 73373, 73148, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 73733, 73508, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 74093, 73868, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 74453, 74228, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 74813, 74588, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 75173, 74948, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 75533, 75308, 15, 1, nmax);

        simdtrf::transform_g_inner(buffer, 75983, 75668, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 76487, 76172, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 76991, 76676, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 77495, 77180, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 77999, 77684, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 78503, 78188, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 79007, 78692, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 79511, 79196, 21, 1, nmax);

        simdtrf::transform_g_inner(buffer, 80120, 79700, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 80792, 80372, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 81464, 81044, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 82136, 81716, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 82808, 82388, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 83480, 83060, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 84152, 83732, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 84824, 84404, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 85616, 85076, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 86480, 85940, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 87344, 86804, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 88208, 87668, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 89072, 88532, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 89936, 89396, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 90800, 90260, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 91664, 91124, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 92663, 91988, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 93743, 93068, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 94823, 94148, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 95903, 95228, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 96983, 96308, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 98063, 97388, 45, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 98468, 73013, 74093, 75983, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 98873, 73373, 74093, 76487, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 99278, 73733, 74093, 76991, 9,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 99683, 74093, 77495, 9, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 100088, 74453, 75533, 77999, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 100493, 74813, 75533, 78503, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 100898, 75173, 75533, 79007, 9,
                                          nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 101303, 75533, 79511, 9, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 101708, 75983, 77495, 80120, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 102275, 76487, 77495, 80792, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 102842, 76991, 77495, 81464, 9,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 103409, 77495, 82136, 9, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 103976, 77999, 79511, 82808, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 104543, 78503, 79511, 83480, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 105110, 79007, 79511, 84152, 9,
                                          nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 105677, 79511, 84824, 9, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 106244, 80120, 82136, 85616, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 107000, 80792, 82136, 86480, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 107756, 81464, 82136, 87344, 9,
                                          nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 108512, 82136, 88208, 9, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 109268, 82808, 84824, 89072, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 110024, 83480, 84824, 89936, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 110780, 84152, 84824, 90800, 9,
                                          nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 111536, 84824, 91664, 9, nmax);

        simdtrf::compute_hrr_geom_010x_pk(buffer, coordinates, 112292, 85616, 88208, 92663, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pk(buffer, coordinates, 113264, 86480, 88208, 93743, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pk(buffer, coordinates, 114236, 87344, 88208, 94823, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_pk(buffer, coordinates, 115208, 89072, 91664, 95903, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_pk(buffer, coordinates, 116180, 89936, 91664, 96983, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_pk(buffer, coordinates, 117152, 90800, 91664, 98063, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 118124, 98468, 99683, 101708, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 118934, 98873, 99683, 102275, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 119744, 99278, 99683, 102842, 9,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 120554, 99683, 103409, 9, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 121364, 100088, 101303, 103976, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 122174, 100493, 101303, 104543, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 122984, 100898, 101303, 105110, 9,
                                          nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 123794, 101303, 105677, 9, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 124604, 101708, 103409, 106244, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 125738, 102275, 103409, 107000, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 126872, 102842, 103409, 107756, 9,
                                          nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 128006, 103409, 108512, 9, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 129140, 103976, 105677, 109268, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 130274, 104543, 105677, 110024, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 131408, 105110, 105677, 110780, 9,
                                          nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 132542, 105677, 111536, 9, nmax);

        simdtrf::compute_hrr_geom_010x_di(buffer, coordinates, 133676, 106244, 108512, 112292, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_di(buffer, coordinates, 135188, 107000, 108512, 113264, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_di(buffer, coordinates, 136700, 107756, 108512, 114236, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_di(buffer, coordinates, 138212, 109268, 111536, 115208, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_di(buffer, coordinates, 139724, 110024, 111536, 116180, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_di(buffer, coordinates, 141236, 110780, 111536, 117152, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 142748, 118124, 120554, 124604, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 144098, 118934, 120554, 125738, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 145448, 119744, 120554, 126872, 9,
                                          nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 146798, 120554, 128006, 9, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 148148, 121364, 123794, 129140, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 149498, 122174, 123794, 130274, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 150848, 122984, 123794, 131408, 9,
                                          nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 152198, 123794, 132542, 9, nmax);

        simdtrf::compute_hrr_geom_010x_fh(buffer, coordinates, 153548, 124604, 128006, 133676, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fh(buffer, coordinates, 155438, 125738, 128006, 135188, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fh(buffer, coordinates, 157328, 126872, 128006, 136700, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fh(buffer, coordinates, 159218, 129140, 132542, 138212, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fh(buffer, coordinates, 161108, 130274, 132542, 139724, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fh(buffer, coordinates, 162998, 131408, 132542, 141236, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gg(buffer, coordinates, 164888, 142748, 146798, 153548, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_gg(buffer, coordinates, 166913, 144098, 146798, 155438, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_gg(buffer, coordinates, 168938, 145448, 146798, 157328, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gg(buffer, coordinates, 170963, 148148, 152198, 159218, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_gg(buffer, coordinates, 172988, 149498, 152198, 161108, 9,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_gg(buffer, coordinates, 175013, 150848, 152198, 162998, 9,
                                          nmax);

        simdtrf::transform_g_inner(buffer, 177038, 170963, 15, 9, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 177038, 81, nmax);

        simdtrf::transform_g_inner(buffer, 177038, 172988, 15, 9, nmax);

        simdtrf::transform_g_outer(values + 729 * nvalues + n * npairs, nvalues, buffer, 177038,
                                   81, nmax);

        simdtrf::transform_g_inner(buffer, 177038, 175013, 15, 9, nmax);

        simdtrf::transform_g_outer(values + 1458 * nvalues + n * npairs, nvalues, buffer, 177038,
                                   81, nmax);

        simdtrf::transform_g_inner(buffer, 177038, 164888, 15, 9, nmax);

        simdtrf::transform_g_outer(values + 2187 * nvalues + n * npairs, nvalues, buffer, 177038,
                                   81, nmax);

        simdtrf::transform_g_inner(buffer, 177038, 166913, 15, 9, nmax);

        simdtrf::transform_g_outer(values + 2916 * nvalues + n * npairs, nvalues, buffer, 177038,
                                   81, nmax);

        simdtrf::transform_g_inner(buffer, 177038, 168938, 15, 9, nmax);

        simdtrf::transform_g_outer(values + 3645 * nvalues + n * npairs, nvalues, buffer, 177038,
                                   81, nmax);
    }

    for (size_t m = 0; m < 4374; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
