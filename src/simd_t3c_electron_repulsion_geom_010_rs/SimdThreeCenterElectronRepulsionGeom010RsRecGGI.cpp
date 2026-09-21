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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGGI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
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
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_ggi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_ggi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 344083, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 6318 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 344083, 184968, 43285, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 15,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 23, 3, 15,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 61, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 64, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 67, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 70, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 73, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 76, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 79, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 82, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 85, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 88, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 91, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 94, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 97, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 100, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 103, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 106, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 109, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 112, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 115, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 118, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 121, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 124, 0, 3, 37, 38,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 127, 0, 3, 38, 39,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 130, 0, 3, 7, 8,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 136, 0, 3, 8, 9,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 142, 0, 3, 9, 10,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 148, 0, 3, 10, 11,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 154, 0, 3, 11, 12,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 160, 0, 3, 12, 13,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 166, 0, 3, 13, 14,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 172, 0, 3, 14, 15,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 178, 0, 3, 15, 16,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 184, 0, 3, 16, 17,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 190, 0, 3, 17, 18,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 196, 0, 3, 18, 19,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 202, 0, 3, 19, 20,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 208, 0, 3, 20, 21,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 214, 0, 3, 24, 25,
                                                                       85, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 220, 0, 3, 25, 26,
                                                                       88, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 226, 0, 3, 26, 27,
                                                                       91, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 232, 0, 3, 27, 28,
                                                                       94, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 238, 0, 3, 28, 29,
                                                                       97, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 244, 0, 3, 29, 30,
                                                                       100, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 250, 0, 3, 30, 31,
                                                                       103, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 256, 0, 3, 31, 32,
                                                                       106, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 262, 0, 3, 32, 33,
                                                                       109, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 268, 0, 3, 33, 34,
                                                                       112, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 274, 0, 3, 34, 35,
                                                                       115, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 280, 0, 3, 35, 36,
                                                                       118, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 286, 0, 3, 36, 37,
                                                                       121, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 292, 0, 3, 37, 38,
                                                                       124, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 40, 43,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 43, 46,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 46, 49,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 49, 52,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 52, 55,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 55, 58,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 58, 61,
                                                                       166, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 61, 64,
                                                                       172, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 64, 67,
                                                                       178, 184, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 67, 70,
                                                                       184, 190, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 70, 73,
                                                                       190, 196, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 73, 76,
                                                                       196, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 76, 79,
                                                                       202, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 85, 88,
                                                                       214, 220, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 88, 91,
                                                                       220, 226, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 91, 94,
                                                                       226, 232, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 94, 97,
                                                                       232, 238, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 97,
                                                                       100, 238, 244, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 478, 0, 3, 100,
                                                                       103, 244, 250, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 488, 0, 3, 103,
                                                                       106, 250, 256, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 106,
                                                                       109, 256, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 508, 0, 3, 109,
                                                                       112, 262, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 112,
                                                                       115, 268, 274, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 528, 0, 3, 115,
                                                                       118, 274, 280, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 538, 0, 3, 118,
                                                                       121, 280, 286, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 548, 0, 3, 121,
                                                                       124, 286, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 558, 0, 3, 130,
                                                                       136, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 573, 0, 3, 136,
                                                                       142, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 588, 0, 3, 142,
                                                                       148, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 603, 0, 3, 148,
                                                                       154, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 618, 0, 3, 154,
                                                                       160, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 633, 0, 3, 160,
                                                                       166, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 648, 0, 3, 166,
                                                                       172, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 663, 0, 3, 172,
                                                                       178, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 678, 0, 3, 178,
                                                                       184, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 693, 0, 3, 184,
                                                                       190, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 708, 0, 3, 190,
                                                                       196, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 723, 0, 3, 196,
                                                                       202, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 738, 0, 3, 214,
                                                                       220, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 753, 0, 3, 220,
                                                                       226, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 768, 0, 3, 226,
                                                                       232, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 783, 0, 3, 232,
                                                                       238, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 798, 0, 3, 238,
                                                                       244, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 813, 0, 3, 244,
                                                                       250, 478, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 828, 0, 3, 250,
                                                                       256, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 843, 0, 3, 256,
                                                                       262, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 858, 0, 3, 262,
                                                                       268, 508, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 873, 0, 3, 268,
                                                                       274, 518, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 888, 0, 3, 274,
                                                                       280, 528, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 903, 0, 3, 280,
                                                                       286, 538, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 918, 0, 3, 298,
                                                                       308, 558, 573, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 939, 0, 3, 308,
                                                                       318, 573, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 960, 0, 3, 318,
                                                                       328, 588, 603, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 981, 0, 3, 328,
                                                                       338, 603, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 338,
                                                                       348, 618, 633, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1023, 0, 3, 348,
                                                                       358, 633, 648, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 358,
                                                                       368, 648, 663, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1065, 0, 3, 368,
                                                                       378, 663, 678, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1086, 0, 3, 378,
                                                                       388, 678, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1107, 0, 3, 388,
                                                                       398, 693, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 398,
                                                                       408, 708, 723, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1149, 0, 3, 428,
                                                                       438, 738, 753, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1170, 0, 3, 438,
                                                                       448, 753, 768, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1191, 0, 3, 448,
                                                                       458, 768, 783, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 458,
                                                                       468, 783, 798, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1233, 0, 3, 468,
                                                                       478, 798, 813, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1254, 0, 3, 478,
                                                                       488, 813, 828, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1275, 0, 3, 488,
                                                                       498, 828, 843, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 498,
                                                                       508, 843, 858, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1317, 0, 3, 508,
                                                                       518, 858, 873, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1338, 0, 3, 518,
                                                                       528, 873, 888, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1359, 0, 3, 528,
                                                                       538, 888, 903, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 558,
                                                                       573, 918, 939, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 573,
                                                                       588, 939, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 588,
                                                                       603, 960, 981, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1464, 0, 3, 603,
                                                                       618, 981, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 618,
                                                                       633, 1002, 1023, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 633,
                                                                       648, 1023, 1044, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 648,
                                                                       663, 1044, 1065, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 663,
                                                                       678, 1065, 1086, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 678,
                                                                       693, 1086, 1107, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 693,
                                                                       708, 1107, 1128, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 738,
                                                                       753, 1149, 1170, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 753,
                                                                       768, 1170, 1191, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 768,
                                                                       783, 1191, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 783,
                                                                       798, 1212, 1233, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 798,
                                                                       813, 1233, 1254, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 813,
                                                                       828, 1254, 1275, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 828,
                                                                       843, 1275, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 843,
                                                                       858, 1296, 1317, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 858,
                                                                       873, 1317, 1338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1912, 0, 3, 873,
                                                                       888, 1338, 1359, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1940, 0, 3, 918,
                                                                       939, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1976, 0, 3, 939,
                                                                       960, 1408, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2012, 0, 3, 960,
                                                                       981, 1436, 1464, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2048, 0, 3, 981,
                                                                       1002, 1464, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2084, 0, 3, 1002,
                                                                       1023, 1492, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2120, 0, 3, 1023,
                                                                       1044, 1520, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2156, 0, 3, 1044,
                                                                       1065, 1548, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2192, 0, 3, 1065,
                                                                       1086, 1576, 1604, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2228, 0, 3, 1086,
                                                                       1107, 1604, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2264, 0, 3, 1149,
                                                                       1170, 1660, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2300, 0, 3, 1170,
                                                                       1191, 1688, 1716, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2336, 0, 3, 1191,
                                                                       1212, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2372, 0, 3, 1212,
                                                                       1233, 1744, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2408, 0, 3, 1233,
                                                                       1254, 1772, 1800, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2444, 0, 3, 1254,
                                                                       1275, 1800, 1828, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2480, 0, 3, 1275,
                                                                       1296, 1828, 1856, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2516, 0, 3, 1296,
                                                                       1317, 1856, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2552, 0, 3, 1317,
                                                                       1338, 1884, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2588, 0, 3, 1380,
                                                                       1408, 1940, 1976, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2633, 0, 3, 1408,
                                                                       1436, 1976, 2012, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2678, 0, 3, 1436,
                                                                       1464, 2012, 2048, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2723, 0, 3, 1464,
                                                                       1492, 2048, 2084, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2768, 0, 3, 1492,
                                                                       1520, 2084, 2120, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2813, 0, 3, 1520,
                                                                       1548, 2120, 2156, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2858, 0, 3, 1548,
                                                                       1576, 2156, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2903, 0, 3, 1576,
                                                                       1604, 2192, 2228, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2948, 0, 3, 1660,
                                                                       1688, 2264, 2300, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2993, 0, 3, 1688,
                                                                       1716, 2300, 2336, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3038, 0, 3, 1716,
                                                                       1744, 2336, 2372, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3083, 0, 3, 1744,
                                                                       1772, 2372, 2408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3128, 0, 3, 1772,
                                                                       1800, 2408, 2444, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3173, 0, 3, 1800,
                                                                       1828, 2444, 2480, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3218, 0, 3, 1828,
                                                                       1856, 2480, 2516, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3263, 0, 3, 1856,
                                                                       1884, 2516, 2552, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3308, 0, 3, 1940,
                                                                       1976, 2588, 2633, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3363, 0, 3, 1976,
                                                                       2012, 2633, 2678, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3418, 0, 3, 2012,
                                                                       2048, 2678, 2723, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3473, 0, 3, 2048,
                                                                       2084, 2723, 2768, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3528, 0, 3, 2084,
                                                                       2120, 2768, 2813, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3583, 0, 3, 2120,
                                                                       2156, 2813, 2858, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3638, 0, 3, 2156,
                                                                       2192, 2858, 2903, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3693, 0, 3, 2264,
                                                                       2300, 2948, 2993, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3748, 0, 3, 2300,
                                                                       2336, 2993, 3038, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3803, 0, 3, 2336,
                                                                       2372, 3038, 3083, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3858, 0, 3, 2372,
                                                                       2408, 3083, 3128, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3913, 0, 3, 2408,
                                                                       2444, 3128, 3173, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 2444,
                                                                       2480, 3173, 3218, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4023, 0, 3, 2480,
                                                                       2516, 3218, 3263, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4078, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4081, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4084, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4087, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4090, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4093, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4096, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4099, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4102, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4105, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4108, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4111, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4114, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4117, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4120, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4123, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4126, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4129, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4132, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4135, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4138, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4141, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4144, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4147, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4150, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4153, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4156, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4159, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4162, 3, 9, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4171, 3, 10, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4180, 3, 11, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4189, 3, 12, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4198, 3, 13, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4207, 3, 14, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4216, 3, 15, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4225, 3, 16, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4234, 3, 17, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4243, 3, 18, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4252, 3, 19, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4261, 3, 20, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4270, 3, 21, 82,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4279, 3, 26, 91,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4288, 3, 27, 94,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4297, 3, 28, 97,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4306, 3, 29, 100,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4315, 3, 30, 103,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4324, 3, 31, 106,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4333, 3, 32, 109,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4342, 3, 33, 112,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4351, 3, 34, 115,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4360, 3, 35, 118,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4369, 3, 36, 121,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4378, 3, 37, 124,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4387, 3, 38, 127,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4396, 3, 46, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4414, 3, 49, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4432, 3, 52, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4450, 3, 55, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4468, 3, 58, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4486, 3, 61, 172,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4504, 3, 64, 178,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4522, 3, 67, 184,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4540, 3, 70, 190,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4558, 3, 73, 196,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4576, 3, 76, 202,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4594, 3, 79, 208,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4612, 3, 91, 226,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4630, 3, 94, 232,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4648, 3, 97, 238,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4666, 3, 100, 244,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4684, 3, 103, 250,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4702, 3, 106, 256,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4720, 3, 109, 262,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4738, 3, 112, 268,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4756, 3, 115, 274,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4774, 3, 118, 280,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4792, 3, 121, 286,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4810, 3, 124, 292,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4828, 3, 142, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4858, 3, 148, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4888, 3, 154, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4918, 3, 160, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4948, 3, 166, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4978, 3, 172, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5008, 3, 178, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5038, 3, 184, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5068, 3, 190, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5098, 3, 196, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5128, 3, 202, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5158, 3, 226, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5188, 3, 232, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5218, 3, 238, 468,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5248, 3, 244, 478,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5278, 3, 250, 488,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5308, 3, 256, 498,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5338, 3, 262, 508,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5368, 3, 268, 518,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5398, 3, 274, 528,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5428, 3, 280, 538,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5458, 3, 286, 548,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5488, 3, 318, 588,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5533, 3, 328, 603,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5578, 3, 338, 618,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5623, 3, 348, 633,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5668, 3, 358, 648,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5713, 3, 368, 663,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5758, 3, 378, 678,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5803, 3, 388, 693,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5848, 3, 398, 708,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5893, 3, 408, 723,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5938, 3, 448, 768,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5983, 3, 458, 783,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6028, 3, 468, 798,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6073, 3, 478, 813,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6118, 3, 488, 828,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6163, 3, 498, 843,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6208, 3, 508, 858,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6253, 3, 518, 873,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6298, 3, 528, 888,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6343, 3, 538, 903,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6388, 3, 588, 960,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6451, 3, 603, 981,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6514, 3, 618,
                                                                       1002, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6577, 3, 633,
                                                                       1023, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6640, 3, 648,
                                                                       1044, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6703, 3, 663,
                                                                       1065, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6766, 3, 678,
                                                                       1086, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6829, 3, 693,
                                                                       1107, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6892, 3, 708,
                                                                       1128, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6955, 3, 768,
                                                                       1191, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7018, 3, 783,
                                                                       1212, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7081, 3, 798,
                                                                       1233, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7144, 3, 813,
                                                                       1254, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7207, 3, 828,
                                                                       1275, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7270, 3, 843,
                                                                       1296, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7333, 3, 858,
                                                                       1317, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7396, 3, 873,
                                                                       1338, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7459, 3, 888,
                                                                       1359, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7522, 3, 960,
                                                                       1436, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7606, 3, 981,
                                                                       1464, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7690, 3, 1002,
                                                                       1492, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7774, 3, 1023,
                                                                       1520, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7858, 3, 1044,
                                                                       1548, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7942, 3, 1065,
                                                                       1576, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8026, 3, 1086,
                                                                       1604, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8110, 3, 1107,
                                                                       1632, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8194, 3, 1191,
                                                                       1716, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8278, 3, 1212,
                                                                       1744, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8362, 3, 1233,
                                                                       1772, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8446, 3, 1254,
                                                                       1800, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8530, 3, 1275,
                                                                       1828, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8614, 3, 1296,
                                                                       1856, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8698, 3, 1317,
                                                                       1884, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8782, 3, 1338,
                                                                       1912, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8866, 3, 1436,
                                                                       2012, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8974, 3, 1464,
                                                                       2048, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9082, 3, 1492,
                                                                       2084, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9190, 3, 1520,
                                                                       2120, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9298, 3, 1548,
                                                                       2156, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9406, 3, 1576,
                                                                       2192, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9514, 3, 1604,
                                                                       2228, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9622, 3, 1716,
                                                                       2336, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9730, 3, 1744,
                                                                       2372, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9838, 3, 1772,
                                                                       2408, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9946, 3, 1800,
                                                                       2444, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10054, 3, 1828,
                                                                       2480, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10162, 3, 1856,
                                                                       2516, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10270, 3, 1884,
                                                                       2552, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10378, 3, 2012,
                                                                       2678, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10513, 3, 2048,
                                                                       2723, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10648, 3, 2084,
                                                                       2768, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10783, 3, 2120,
                                                                       2813, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10918, 3, 2156,
                                                                       2858, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11053, 3, 2192,
                                                                       2903, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11188, 3, 2336,
                                                                       3038, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11323, 3, 2372,
                                                                       3083, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11458, 3, 2408,
                                                                       3128, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11593, 3, 2444,
                                                                       3173, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11728, 3, 2480,
                                                                       3218, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11863, 3, 2516,
                                                                       3263, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11998, 3, 2678,
                                                                       3418, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12163, 3, 2723,
                                                                       3473, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12328, 3, 2768,
                                                                       3528, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12493, 3, 2813,
                                                                       3583, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12658, 3, 2858,
                                                                       3638, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12823, 3, 3038,
                                                                       3803, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12988, 3, 3083,
                                                                       3858, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13153, 3, 3128,
                                                                       3913, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13318, 3, 3173,
                                                                       3968, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13483, 3, 3218,
                                                                       4023, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13648, 3, 7, 8,
                                                                       4078, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13654, 3, 8, 9,
                                                                       4081, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13660, 3, 9, 10,
                                                                       4084, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13666, 3, 10, 11,
                                                                       4087, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13672, 3, 11, 12,
                                                                       4090, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13678, 3, 12, 13,
                                                                       4093, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13684, 3, 13, 14,
                                                                       4096, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13690, 3, 14, 15,
                                                                       4099, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13696, 3, 15, 16,
                                                                       4102, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13702, 3, 16, 17,
                                                                       4105, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13708, 3, 17, 18,
                                                                       4108, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13714, 3, 18, 19,
                                                                       4111, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13720, 3, 19, 20,
                                                                       4114, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13726, 3, 20, 21,
                                                                       4117, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13732, 3, 24, 25,
                                                                       4120, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13738, 3, 25, 26,
                                                                       4123, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13744, 3, 26, 27,
                                                                       4126, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13750, 3, 27, 28,
                                                                       4129, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13756, 3, 28, 29,
                                                                       4132, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13762, 3, 29, 30,
                                                                       4135, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13768, 3, 30, 31,
                                                                       4138, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13774, 3, 31, 32,
                                                                       4141, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13780, 3, 32, 33,
                                                                       4144, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13786, 3, 33, 34,
                                                                       4147, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13792, 3, 34, 35,
                                                                       4150, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13798, 3, 35, 36,
                                                                       4153, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13804, 3, 36, 37,
                                                                       4156, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 13810, 3, 37, 38,
                                                                       4159, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13816, 0, 3,
                                                                       13648, 4078, 13654, 4162,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13834, 0, 3,
                                                                       13654, 4081, 13660, 4171,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13852, 0, 3,
                                                                       13660, 4084, 13666, 4180,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13870, 0, 3,
                                                                       13666, 4087, 13672, 4189,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13888, 0, 3,
                                                                       13672, 4090, 13678, 4198,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13906, 0, 3,
                                                                       13678, 4093, 13684, 4207,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13924, 0, 3,
                                                                       13684, 4096, 13690, 4216,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13942, 0, 3,
                                                                       13690, 4099, 13696, 4225,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13960, 0, 3,
                                                                       13696, 4102, 13702, 4234,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13978, 0, 3,
                                                                       13702, 4105, 13708, 4243,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13996, 0, 3,
                                                                       13708, 4108, 13714, 4252,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14014, 0, 3,
                                                                       13714, 4111, 13720, 4261,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14032, 0, 3,
                                                                       13720, 4114, 13726, 4270,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14050, 0, 3,
                                                                       13732, 4120, 13738, 4279,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14068, 0, 3,
                                                                       13738, 4123, 13744, 4288,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14086, 0, 3,
                                                                       13744, 4126, 13750, 4297,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14104, 0, 3,
                                                                       13750, 4129, 13756, 4306,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14122, 0, 3,
                                                                       13756, 4132, 13762, 4315,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14140, 0, 3,
                                                                       13762, 4135, 13768, 4324,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14158, 0, 3,
                                                                       13768, 4138, 13774, 4333,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14176, 0, 3,
                                                                       13774, 4141, 13780, 4342,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14194, 0, 3,
                                                                       13780, 4144, 13786, 4351,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14212, 0, 3,
                                                                       13786, 4147, 13792, 4360,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14230, 0, 3,
                                                                       13792, 4150, 13798, 4369,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14248, 0, 3,
                                                                       13798, 4153, 13804, 4378,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14266, 0, 3,
                                                                       13804, 4156, 13810, 4387,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14284, 0, 3,
                                                                       13816, 4162, 13834, 130,
                                                                       136, 4396, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14320, 0, 3,
                                                                       13834, 4171, 13852, 136,
                                                                       142, 4414, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14356, 0, 3,
                                                                       13852, 4180, 13870, 142,
                                                                       148, 4432, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14392, 0, 3,
                                                                       13870, 4189, 13888, 148,
                                                                       154, 4450, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14428, 0, 3,
                                                                       13888, 4198, 13906, 154,
                                                                       160, 4468, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14464, 0, 3,
                                                                       13906, 4207, 13924, 160,
                                                                       166, 4486, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14500, 0, 3,
                                                                       13924, 4216, 13942, 166,
                                                                       172, 4504, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14536, 0, 3,
                                                                       13942, 4225, 13960, 172,
                                                                       178, 4522, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14572, 0, 3,
                                                                       13960, 4234, 13978, 178,
                                                                       184, 4540, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14608, 0, 3,
                                                                       13978, 4243, 13996, 184,
                                                                       190, 4558, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14644, 0, 3,
                                                                       13996, 4252, 14014, 190,
                                                                       196, 4576, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14680, 0, 3,
                                                                       14014, 4261, 14032, 196,
                                                                       202, 4594, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14716, 0, 3,
                                                                       14050, 4279, 14068, 214,
                                                                       220, 4612, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14752, 0, 3,
                                                                       14068, 4288, 14086, 220,
                                                                       226, 4630, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14788, 0, 3,
                                                                       14086, 4297, 14104, 226,
                                                                       232, 4648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14824, 0, 3,
                                                                       14104, 4306, 14122, 232,
                                                                       238, 4666, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14860, 0, 3,
                                                                       14122, 4315, 14140, 238,
                                                                       244, 4684, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14896, 0, 3,
                                                                       14140, 4324, 14158, 244,
                                                                       250, 4702, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14932, 0, 3,
                                                                       14158, 4333, 14176, 250,
                                                                       256, 4720, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14968, 0, 3,
                                                                       14176, 4342, 14194, 256,
                                                                       262, 4738, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15004, 0, 3,
                                                                       14194, 4351, 14212, 262,
                                                                       268, 4756, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15040, 0, 3,
                                                                       14212, 4360, 14230, 268,
                                                                       274, 4774, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15076, 0, 3,
                                                                       14230, 4369, 14248, 274,
                                                                       280, 4792, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15112, 0, 3,
                                                                       14248, 4378, 14266, 280,
                                                                       286, 4810, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15148, 0, 3,
                                                                       14284, 4396, 14320, 298,
                                                                       308, 4828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15208, 0, 3,
                                                                       14320, 4414, 14356, 308,
                                                                       318, 4858, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15268, 0, 3,
                                                                       14356, 4432, 14392, 318,
                                                                       328, 4888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15328, 0, 3,
                                                                       14392, 4450, 14428, 328,
                                                                       338, 4918, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15388, 0, 3,
                                                                       14428, 4468, 14464, 338,
                                                                       348, 4948, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15448, 0, 3,
                                                                       14464, 4486, 14500, 348,
                                                                       358, 4978, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15508, 0, 3,
                                                                       14500, 4504, 14536, 358,
                                                                       368, 5008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15568, 0, 3,
                                                                       14536, 4522, 14572, 368,
                                                                       378, 5038, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15628, 0, 3,
                                                                       14572, 4540, 14608, 378,
                                                                       388, 5068, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15688, 0, 3,
                                                                       14608, 4558, 14644, 388,
                                                                       398, 5098, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15748, 0, 3,
                                                                       14644, 4576, 14680, 398,
                                                                       408, 5128, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15808, 0, 3,
                                                                       14716, 4612, 14752, 428,
                                                                       438, 5158, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15868, 0, 3,
                                                                       14752, 4630, 14788, 438,
                                                                       448, 5188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15928, 0, 3,
                                                                       14788, 4648, 14824, 448,
                                                                       458, 5218, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15988, 0, 3,
                                                                       14824, 4666, 14860, 458,
                                                                       468, 5248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16048, 0, 3,
                                                                       14860, 4684, 14896, 468,
                                                                       478, 5278, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16108, 0, 3,
                                                                       14896, 4702, 14932, 478,
                                                                       488, 5308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16168, 0, 3,
                                                                       14932, 4720, 14968, 488,
                                                                       498, 5338, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16228, 0, 3,
                                                                       14968, 4738, 15004, 498,
                                                                       508, 5368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16288, 0, 3,
                                                                       15004, 4756, 15040, 508,
                                                                       518, 5398, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16348, 0, 3,
                                                                       15040, 4774, 15076, 518,
                                                                       528, 5428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16408, 0, 3,
                                                                       15076, 4792, 15112, 528,
                                                                       538, 5458, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16468, 0, 3,
                                                                       15148, 4828, 15208, 558,
                                                                       573, 5488, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16558, 0, 3,
                                                                       15208, 4858, 15268, 573,
                                                                       588, 5533, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16648, 0, 3,
                                                                       15268, 4888, 15328, 588,
                                                                       603, 5578, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16738, 0, 3,
                                                                       15328, 4918, 15388, 603,
                                                                       618, 5623, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16828, 0, 3,
                                                                       15388, 4948, 15448, 618,
                                                                       633, 5668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16918, 0, 3,
                                                                       15448, 4978, 15508, 633,
                                                                       648, 5713, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17008, 0, 3,
                                                                       15508, 5008, 15568, 648,
                                                                       663, 5758, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17098, 0, 3,
                                                                       15568, 5038, 15628, 663,
                                                                       678, 5803, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17188, 0, 3,
                                                                       15628, 5068, 15688, 678,
                                                                       693, 5848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17278, 0, 3,
                                                                       15688, 5098, 15748, 693,
                                                                       708, 5893, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17368, 0, 3,
                                                                       15808, 5158, 15868, 738,
                                                                       753, 5938, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17458, 0, 3,
                                                                       15868, 5188, 15928, 753,
                                                                       768, 5983, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17548, 0, 3,
                                                                       15928, 5218, 15988, 768,
                                                                       783, 6028, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17638, 0, 3,
                                                                       15988, 5248, 16048, 783,
                                                                       798, 6073, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17728, 0, 3,
                                                                       16048, 5278, 16108, 798,
                                                                       813, 6118, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17818, 0, 3,
                                                                       16108, 5308, 16168, 813,
                                                                       828, 6163, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17908, 0, 3,
                                                                       16168, 5338, 16228, 828,
                                                                       843, 6208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17998, 0, 3,
                                                                       16228, 5368, 16288, 843,
                                                                       858, 6253, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 18088, 0, 3,
                                                                       16288, 5398, 16348, 858,
                                                                       873, 6298, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 18178, 0, 3,
                                                                       16348, 5428, 16408, 873,
                                                                       888, 6343, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18268, 0, 3,
                                                                       16468, 5488, 16558, 918,
                                                                       939, 6388, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18394, 0, 3,
                                                                       16558, 5533, 16648, 939,
                                                                       960, 6451, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18520, 0, 3,
                                                                       16648, 5578, 16738, 960,
                                                                       981, 6514, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18646, 0, 3,
                                                                       16738, 5623, 16828, 981,
                                                                       1002, 6577, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18772, 0, 3,
                                                                       16828, 5668, 16918, 1002,
                                                                       1023, 6640, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18898, 0, 3,
                                                                       16918, 5713, 17008, 1023,
                                                                       1044, 6703, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19024, 0, 3,
                                                                       17008, 5758, 17098, 1044,
                                                                       1065, 6766, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19150, 0, 3,
                                                                       17098, 5803, 17188, 1065,
                                                                       1086, 6829, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19276, 0, 3,
                                                                       17188, 5848, 17278, 1086,
                                                                       1107, 6892, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19402, 0, 3,
                                                                       17368, 5938, 17458, 1149,
                                                                       1170, 6955, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19528, 0, 3,
                                                                       17458, 5983, 17548, 1170,
                                                                       1191, 7018, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19654, 0, 3,
                                                                       17548, 6028, 17638, 1191,
                                                                       1212, 7081, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19780, 0, 3,
                                                                       17638, 6073, 17728, 1212,
                                                                       1233, 7144, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19906, 0, 3,
                                                                       17728, 6118, 17818, 1233,
                                                                       1254, 7207, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 20032, 0, 3,
                                                                       17818, 6163, 17908, 1254,
                                                                       1275, 7270, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 20158, 0, 3,
                                                                       17908, 6208, 17998, 1275,
                                                                       1296, 7333, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 20284, 0, 3,
                                                                       17998, 6253, 18088, 1296,
                                                                       1317, 7396, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 20410, 0, 3,
                                                                       18088, 6298, 18178, 1317,
                                                                       1338, 7459, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20536, 0, 3,
                                                                       18268, 6388, 18394, 1380,
                                                                       1408, 7522, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20704, 0, 3,
                                                                       18394, 6451, 18520, 1408,
                                                                       1436, 7606, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20872, 0, 3,
                                                                       18520, 6514, 18646, 1436,
                                                                       1464, 7690, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21040, 0, 3,
                                                                       18646, 6577, 18772, 1464,
                                                                       1492, 7774, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21208, 0, 3,
                                                                       18772, 6640, 18898, 1492,
                                                                       1520, 7858, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21376, 0, 3,
                                                                       18898, 6703, 19024, 1520,
                                                                       1548, 7942, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21544, 0, 3,
                                                                       19024, 6766, 19150, 1548,
                                                                       1576, 8026, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21712, 0, 3,
                                                                       19150, 6829, 19276, 1576,
                                                                       1604, 8110, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21880, 0, 3,
                                                                       19402, 6955, 19528, 1660,
                                                                       1688, 8194, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 22048, 0, 3,
                                                                       19528, 7018, 19654, 1688,
                                                                       1716, 8278, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 22216, 0, 3,
                                                                       19654, 7081, 19780, 1716,
                                                                       1744, 8362, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 22384, 0, 3,
                                                                       19780, 7144, 19906, 1744,
                                                                       1772, 8446, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 22552, 0, 3,
                                                                       19906, 7207, 20032, 1772,
                                                                       1800, 8530, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 22720, 0, 3,
                                                                       20032, 7270, 20158, 1800,
                                                                       1828, 8614, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 22888, 0, 3,
                                                                       20158, 7333, 20284, 1828,
                                                                       1856, 8698, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23056, 0, 3,
                                                                       20284, 7396, 20410, 1856,
                                                                       1884, 8782, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23224, 0, 3,
                                                                       20536, 7522, 20704, 1940,
                                                                       1976, 8866, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23440, 0, 3,
                                                                       20704, 7606, 20872, 1976,
                                                                       2012, 8974, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23656, 0, 3,
                                                                       20872, 7690, 21040, 2012,
                                                                       2048, 9082, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23872, 0, 3,
                                                                       21040, 7774, 21208, 2048,
                                                                       2084, 9190, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 24088, 0, 3,
                                                                       21208, 7858, 21376, 2084,
                                                                       2120, 9298, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 24304, 0, 3,
                                                                       21376, 7942, 21544, 2120,
                                                                       2156, 9406, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 24520, 0, 3,
                                                                       21544, 8026, 21712, 2156,
                                                                       2192, 9514, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 24736, 0, 3,
                                                                       21880, 8194, 22048, 2264,
                                                                       2300, 9622, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 24952, 0, 3,
                                                                       22048, 8278, 22216, 2300,
                                                                       2336, 9730, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 25168, 0, 3,
                                                                       22216, 8362, 22384, 2336,
                                                                       2372, 9838, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 25384, 0, 3,
                                                                       22384, 8446, 22552, 2372,
                                                                       2408, 9946, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 25600, 0, 3,
                                                                       22552, 8530, 22720, 2408,
                                                                       2444, 10054, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 25816, 0, 3,
                                                                       22720, 8614, 22888, 2444,
                                                                       2480, 10162, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 26032, 0, 3,
                                                                       22888, 8698, 23056, 2480,
                                                                       2516, 10270, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 26248, 0, 3,
                                                                       23224, 8866, 23440, 2588,
                                                                       2633, 10378, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 26518, 0, 3,
                                                                       23440, 8974, 23656, 2633,
                                                                       2678, 10513, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 26788, 0, 3,
                                                                       23656, 9082, 23872, 2678,
                                                                       2723, 10648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 27058, 0, 3,
                                                                       23872, 9190, 24088, 2723,
                                                                       2768, 10783, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 27328, 0, 3,
                                                                       24088, 9298, 24304, 2768,
                                                                       2813, 10918, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 27598, 0, 3,
                                                                       24304, 9406, 24520, 2813,
                                                                       2858, 11053, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 27868, 0, 3,
                                                                       24736, 9622, 24952, 2948,
                                                                       2993, 11188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 28138, 0, 3,
                                                                       24952, 9730, 25168, 2993,
                                                                       3038, 11323, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 28408, 0, 3,
                                                                       25168, 9838, 25384, 3038,
                                                                       3083, 11458, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 28678, 0, 3,
                                                                       25384, 9946, 25600, 3083,
                                                                       3128, 11593, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 28948, 0, 3,
                                                                       25600, 10054, 25816, 3128,
                                                                       3173, 11728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 29218, 0, 3,
                                                                       25816, 10162, 26032, 3173,
                                                                       3218, 11863, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 29488, 0, 3,
                                                                       26248, 10378, 26518, 3308,
                                                                       3363, 11998, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 29818, 0, 3,
                                                                       26518, 10513, 26788, 3363,
                                                                       3418, 12163, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 30148, 0, 3,
                                                                       26788, 10648, 27058, 3418,
                                                                       3473, 12328, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 30478, 0, 3,
                                                                       27058, 10783, 27328, 3473,
                                                                       3528, 12493, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 30808, 0, 3,
                                                                       27328, 10918, 27598, 3528,
                                                                       3583, 12658, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 31138, 0, 3,
                                                                       27868, 11188, 28138, 3693,
                                                                       3748, 12823, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 31468, 0, 3,
                                                                       28138, 11323, 28408, 3748,
                                                                       3803, 12988, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 31798, 0, 3,
                                                                       28408, 11458, 28678, 3803,
                                                                       3858, 13153, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 32128, 0, 3,
                                                                       28678, 11593, 28948, 3858,
                                                                       3913, 13318, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 32458, 0, 3,
                                                                       28948, 11728, 29218, 3913,
                                                                       3968, 13483, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32788, 3, 4078,
                                                                       4081, 13660, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32798, 3, 4081,
                                                                       4084, 13666, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32808, 3, 4084,
                                                                       4087, 13672, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32818, 3, 4087,
                                                                       4090, 13678, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32828, 3, 4090,
                                                                       4093, 13684, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32838, 3, 4093,
                                                                       4096, 13690, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32848, 3, 4096,
                                                                       4099, 13696, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32858, 3, 4099,
                                                                       4102, 13702, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32868, 3, 4102,
                                                                       4105, 13708, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32878, 3, 4105,
                                                                       4108, 13714, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32888, 3, 4108,
                                                                       4111, 13720, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32898, 3, 4111,
                                                                       4114, 13726, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32908, 3, 4120,
                                                                       4123, 13744, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32918, 3, 4123,
                                                                       4126, 13750, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32928, 3, 4126,
                                                                       4129, 13756, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32938, 3, 4129,
                                                                       4132, 13762, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32948, 3, 4132,
                                                                       4135, 13768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32958, 3, 4135,
                                                                       4138, 13774, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32968, 3, 4138,
                                                                       4141, 13780, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32978, 3, 4141,
                                                                       4144, 13786, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32988, 3, 4144,
                                                                       4147, 13792, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 32998, 3, 4147,
                                                                       4150, 13798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33008, 3, 4150,
                                                                       4153, 13804, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33018, 3, 4153,
                                                                       4156, 13810, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33028, 0, 3,
                                                                       32788, 13660, 32798,
                                                                       13852, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33058, 0, 3,
                                                                       32798, 13666, 32808,
                                                                       13870, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33088, 0, 3,
                                                                       32808, 13672, 32818,
                                                                       13888, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33118, 0, 3,
                                                                       32818, 13678, 32828,
                                                                       13906, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33148, 0, 3,
                                                                       32828, 13684, 32838,
                                                                       13924, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33178, 0, 3,
                                                                       32838, 13690, 32848,
                                                                       13942, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33208, 0, 3,
                                                                       32848, 13696, 32858,
                                                                       13960, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33238, 0, 3,
                                                                       32858, 13702, 32868,
                                                                       13978, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33268, 0, 3,
                                                                       32868, 13708, 32878,
                                                                       13996, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33298, 0, 3,
                                                                       32878, 13714, 32888,
                                                                       14014, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33328, 0, 3,
                                                                       32888, 13720, 32898,
                                                                       14032, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33358, 0, 3,
                                                                       32908, 13744, 32918,
                                                                       14086, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33388, 0, 3,
                                                                       32918, 13750, 32928,
                                                                       14104, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33418, 0, 3,
                                                                       32928, 13756, 32938,
                                                                       14122, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33448, 0, 3,
                                                                       32938, 13762, 32948,
                                                                       14140, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33478, 0, 3,
                                                                       32948, 13768, 32958,
                                                                       14158, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33508, 0, 3,
                                                                       32958, 13774, 32968,
                                                                       14176, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33538, 0, 3,
                                                                       32968, 13780, 32978,
                                                                       14194, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33568, 0, 3,
                                                                       32978, 13786, 32988,
                                                                       14212, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33598, 0, 3,
                                                                       32988, 13792, 32998,
                                                                       14230, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33628, 0, 3,
                                                                       32998, 13798, 33008,
                                                                       14248, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33658, 0, 3,
                                                                       33008, 13804, 33018,
                                                                       14266, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 33688, 0, 3,
                                                                       33028, 13852, 33058, 4396,
                                                                       4414, 14356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 33748, 0, 3,
                                                                       33058, 13870, 33088, 4414,
                                                                       4432, 14392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 33808, 0, 3,
                                                                       33088, 13888, 33118, 4432,
                                                                       4450, 14428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 33868, 0, 3,
                                                                       33118, 13906, 33148, 4450,
                                                                       4468, 14464, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 33928, 0, 3,
                                                                       33148, 13924, 33178, 4468,
                                                                       4486, 14500, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 33988, 0, 3,
                                                                       33178, 13942, 33208, 4486,
                                                                       4504, 14536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34048, 0, 3,
                                                                       33208, 13960, 33238, 4504,
                                                                       4522, 14572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34108, 0, 3,
                                                                       33238, 13978, 33268, 4522,
                                                                       4540, 14608, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34168, 0, 3,
                                                                       33268, 13996, 33298, 4540,
                                                                       4558, 14644, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34228, 0, 3,
                                                                       33298, 14014, 33328, 4558,
                                                                       4576, 14680, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34288, 0, 3,
                                                                       33358, 14086, 33388, 4612,
                                                                       4630, 14788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34348, 0, 3,
                                                                       33388, 14104, 33418, 4630,
                                                                       4648, 14824, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34408, 0, 3,
                                                                       33418, 14122, 33448, 4648,
                                                                       4666, 14860, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34468, 0, 3,
                                                                       33448, 14140, 33478, 4666,
                                                                       4684, 14896, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34528, 0, 3,
                                                                       33478, 14158, 33508, 4684,
                                                                       4702, 14932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34588, 0, 3,
                                                                       33508, 14176, 33538, 4702,
                                                                       4720, 14968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34648, 0, 3,
                                                                       33538, 14194, 33568, 4720,
                                                                       4738, 15004, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34708, 0, 3,
                                                                       33568, 14212, 33598, 4738,
                                                                       4756, 15040, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34768, 0, 3,
                                                                       33598, 14230, 33628, 4756,
                                                                       4774, 15076, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34828, 0, 3,
                                                                       33628, 14248, 33658, 4774,
                                                                       4792, 15112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 34888, 0, 3,
                                                                       33688, 14356, 33748, 4828,
                                                                       4858, 15268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 34988, 0, 3,
                                                                       33748, 14392, 33808, 4858,
                                                                       4888, 15328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35088, 0, 3,
                                                                       33808, 14428, 33868, 4888,
                                                                       4918, 15388, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35188, 0, 3,
                                                                       33868, 14464, 33928, 4918,
                                                                       4948, 15448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35288, 0, 3,
                                                                       33928, 14500, 33988, 4948,
                                                                       4978, 15508, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35388, 0, 3,
                                                                       33988, 14536, 34048, 4978,
                                                                       5008, 15568, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35488, 0, 3,
                                                                       34048, 14572, 34108, 5008,
                                                                       5038, 15628, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35588, 0, 3,
                                                                       34108, 14608, 34168, 5038,
                                                                       5068, 15688, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35688, 0, 3,
                                                                       34168, 14644, 34228, 5068,
                                                                       5098, 15748, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35788, 0, 3,
                                                                       34288, 14788, 34348, 5158,
                                                                       5188, 15928, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35888, 0, 3,
                                                                       34348, 14824, 34408, 5188,
                                                                       5218, 15988, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35988, 0, 3,
                                                                       34408, 14860, 34468, 5218,
                                                                       5248, 16048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36088, 0, 3,
                                                                       34468, 14896, 34528, 5248,
                                                                       5278, 16108, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36188, 0, 3,
                                                                       34528, 14932, 34588, 5278,
                                                                       5308, 16168, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36288, 0, 3,
                                                                       34588, 14968, 34648, 5308,
                                                                       5338, 16228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36388, 0, 3,
                                                                       34648, 15004, 34708, 5338,
                                                                       5368, 16288, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36488, 0, 3,
                                                                       34708, 15040, 34768, 5368,
                                                                       5398, 16348, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36588, 0, 3,
                                                                       34768, 15076, 34828, 5398,
                                                                       5428, 16408, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 36688, 0, 3,
                                                                       34888, 15268, 34988, 5488,
                                                                       5533, 16648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 36838, 0, 3,
                                                                       34988, 15328, 35088, 5533,
                                                                       5578, 16738, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 36988, 0, 3,
                                                                       35088, 15388, 35188, 5578,
                                                                       5623, 16828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 37138, 0, 3,
                                                                       35188, 15448, 35288, 5623,
                                                                       5668, 16918, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 37288, 0, 3,
                                                                       35288, 15508, 35388, 5668,
                                                                       5713, 17008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 37438, 0, 3,
                                                                       35388, 15568, 35488, 5713,
                                                                       5758, 17098, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 37588, 0, 3,
                                                                       35488, 15628, 35588, 5758,
                                                                       5803, 17188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 37738, 0, 3,
                                                                       35588, 15688, 35688, 5803,
                                                                       5848, 17278, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 37888, 0, 3,
                                                                       35788, 15928, 35888, 5938,
                                                                       5983, 17548, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38038, 0, 3,
                                                                       35888, 15988, 35988, 5983,
                                                                       6028, 17638, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38188, 0, 3,
                                                                       35988, 16048, 36088, 6028,
                                                                       6073, 17728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38338, 0, 3,
                                                                       36088, 16108, 36188, 6073,
                                                                       6118, 17818, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38488, 0, 3,
                                                                       36188, 16168, 36288, 6118,
                                                                       6163, 17908, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38638, 0, 3,
                                                                       36288, 16228, 36388, 6163,
                                                                       6208, 17998, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38788, 0, 3,
                                                                       36388, 16288, 36488, 6208,
                                                                       6253, 18088, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38938, 0, 3,
                                                                       36488, 16348, 36588, 6253,
                                                                       6298, 18178, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 39088, 0, 3,
                                                                       36688, 16648, 36838, 6388,
                                                                       6451, 18520, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 39298, 0, 3,
                                                                       36838, 16738, 36988, 6451,
                                                                       6514, 18646, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 39508, 0, 3,
                                                                       36988, 16828, 37138, 6514,
                                                                       6577, 18772, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 39718, 0, 3,
                                                                       37138, 16918, 37288, 6577,
                                                                       6640, 18898, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 39928, 0, 3,
                                                                       37288, 17008, 37438, 6640,
                                                                       6703, 19024, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40138, 0, 3,
                                                                       37438, 17098, 37588, 6703,
                                                                       6766, 19150, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40348, 0, 3,
                                                                       37588, 17188, 37738, 6766,
                                                                       6829, 19276, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40558, 0, 3,
                                                                       37888, 17548, 38038, 6955,
                                                                       7018, 19654, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40768, 0, 3,
                                                                       38038, 17638, 38188, 7018,
                                                                       7081, 19780, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40978, 0, 3,
                                                                       38188, 17728, 38338, 7081,
                                                                       7144, 19906, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41188, 0, 3,
                                                                       38338, 17818, 38488, 7144,
                                                                       7207, 20032, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41398, 0, 3,
                                                                       38488, 17908, 38638, 7207,
                                                                       7270, 20158, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41608, 0, 3,
                                                                       38638, 17998, 38788, 7270,
                                                                       7333, 20284, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41818, 0, 3,
                                                                       38788, 18088, 38938, 7333,
                                                                       7396, 20410, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 42028, 0, 3,
                                                                       39088, 18520, 39298, 7522,
                                                                       7606, 20872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 42308, 0, 3,
                                                                       39298, 18646, 39508, 7606,
                                                                       7690, 21040, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 42588, 0, 3,
                                                                       39508, 18772, 39718, 7690,
                                                                       7774, 21208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 42868, 0, 3,
                                                                       39718, 18898, 39928, 7774,
                                                                       7858, 21376, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43148, 0, 3,
                                                                       39928, 19024, 40138, 7858,
                                                                       7942, 21544, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43428, 0, 3,
                                                                       40138, 19150, 40348, 7942,
                                                                       8026, 21712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43708, 0, 3,
                                                                       40558, 19654, 40768, 8194,
                                                                       8278, 22216, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43988, 0, 3,
                                                                       40768, 19780, 40978, 8278,
                                                                       8362, 22384, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 44268, 0, 3,
                                                                       40978, 19906, 41188, 8362,
                                                                       8446, 22552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 44548, 0, 3,
                                                                       41188, 20032, 41398, 8446,
                                                                       8530, 22720, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 44828, 0, 3,
                                                                       41398, 20158, 41608, 8530,
                                                                       8614, 22888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 45108, 0, 3,
                                                                       41608, 20284, 41818, 8614,
                                                                       8698, 23056, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 45388, 0, 3,
                                                                       42028, 20872, 42308, 8866,
                                                                       8974, 23656, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 45748, 0, 3,
                                                                       42308, 21040, 42588, 8974,
                                                                       9082, 23872, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 46108, 0, 3,
                                                                       42588, 21208, 42868, 9082,
                                                                       9190, 24088, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 46468, 0, 3,
                                                                       42868, 21376, 43148, 9190,
                                                                       9298, 24304, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 46828, 0, 3,
                                                                       43148, 21544, 43428, 9298,
                                                                       9406, 24520, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 47188, 0, 3,
                                                                       43708, 22216, 43988, 9622,
                                                                       9730, 25168, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 47548, 0, 3,
                                                                       43988, 22384, 44268, 9730,
                                                                       9838, 25384, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 47908, 0, 3,
                                                                       44268, 22552, 44548, 9838,
                                                                       9946, 25600, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 48268, 0, 3,
                                                                       44548, 22720, 44828, 9946,
                                                                       10054, 25816, ncols,
                                                                       gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 48628, 0, 3,
                                                                       44828, 22888, 45108,
                                                                       10054, 10162, 26032,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 48988, 0, 3,
                                                                       45388, 23656, 45748,
                                                                       10378, 10513, 26788,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 49438, 0, 3,
                                                                       45748, 23872, 46108,
                                                                       10513, 10648, 27058,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 49888, 0, 3,
                                                                       46108, 24088, 46468,
                                                                       10648, 10783, 27328,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 50338, 0, 3,
                                                                       46468, 24304, 46828,
                                                                       10783, 10918, 27598,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 50788, 0, 3,
                                                                       47188, 25168, 47548,
                                                                       11188, 11323, 28408,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 51238, 0, 3,
                                                                       47548, 25384, 47908,
                                                                       11323, 11458, 28678,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 51688, 0, 3,
                                                                       47908, 25600, 48268,
                                                                       11458, 11593, 28948,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 52138, 0, 3,
                                                                       48268, 25816, 48628,
                                                                       11593, 11728, 29218,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 52588, 0, 3,
                                                                       48988, 26788, 49438,
                                                                       11998, 12163, 30148,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 53138, 0, 3,
                                                                       49438, 27058, 49888,
                                                                       12163, 12328, 30478,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 53688, 0, 3,
                                                                       49888, 27328, 50338,
                                                                       12328, 12493, 30808,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 54238, 0, 3,
                                                                       50788, 28408, 51238,
                                                                       12823, 12988, 31798,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 54788, 0, 3,
                                                                       51238, 28678, 51688,
                                                                       12988, 13153, 32128,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 55338, 0, 3,
                                                                       51688, 28948, 52138,
                                                                       13153, 13318, 32458,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 55888, 3, 13648,
                                                                       13654, 32788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 55903, 3, 13654,
                                                                       13660, 32798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 55918, 3, 13660,
                                                                       13666, 32808, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 55933, 3, 13666,
                                                                       13672, 32818, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 55948, 3, 13672,
                                                                       13678, 32828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 55963, 3, 13678,
                                                                       13684, 32838, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 55978, 3, 13684,
                                                                       13690, 32848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 55993, 3, 13690,
                                                                       13696, 32858, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56008, 3, 13696,
                                                                       13702, 32868, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56023, 3, 13702,
                                                                       13708, 32878, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56038, 3, 13708,
                                                                       13714, 32888, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56053, 3, 13714,
                                                                       13720, 32898, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56068, 3, 13732,
                                                                       13738, 32908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56083, 3, 13738,
                                                                       13744, 32918, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56098, 3, 13744,
                                                                       13750, 32928, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56113, 3, 13750,
                                                                       13756, 32938, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56128, 3, 13756,
                                                                       13762, 32948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56143, 3, 13762,
                                                                       13768, 32958, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56158, 3, 13768,
                                                                       13774, 32968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56173, 3, 13774,
                                                                       13780, 32978, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56188, 3, 13780,
                                                                       13786, 32988, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56203, 3, 13786,
                                                                       13792, 32998, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56218, 3, 13792,
                                                                       13798, 33008, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 56233, 3, 13798,
                                                                       13804, 33018, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56248, 0, 3,
                                                                       55888, 32788, 55903,
                                                                       13816, 13834, 33028,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56293, 0, 3,
                                                                       55903, 32798, 55918,
                                                                       13834, 13852, 33058,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56338, 0, 3,
                                                                       55918, 32808, 55933,
                                                                       13852, 13870, 33088,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56383, 0, 3,
                                                                       55933, 32818, 55948,
                                                                       13870, 13888, 33118,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56428, 0, 3,
                                                                       55948, 32828, 55963,
                                                                       13888, 13906, 33148,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56473, 0, 3,
                                                                       55963, 32838, 55978,
                                                                       13906, 13924, 33178,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56518, 0, 3,
                                                                       55978, 32848, 55993,
                                                                       13924, 13942, 33208,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56563, 0, 3,
                                                                       55993, 32858, 56008,
                                                                       13942, 13960, 33238,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56608, 0, 3,
                                                                       56008, 32868, 56023,
                                                                       13960, 13978, 33268,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56653, 0, 3,
                                                                       56023, 32878, 56038,
                                                                       13978, 13996, 33298,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56698, 0, 3,
                                                                       56038, 32888, 56053,
                                                                       13996, 14014, 33328,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56743, 0, 3,
                                                                       56068, 32908, 56083,
                                                                       14050, 14068, 33358,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56788, 0, 3,
                                                                       56083, 32918, 56098,
                                                                       14068, 14086, 33388,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56833, 0, 3,
                                                                       56098, 32928, 56113,
                                                                       14086, 14104, 33418,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56878, 0, 3,
                                                                       56113, 32938, 56128,
                                                                       14104, 14122, 33448,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56923, 0, 3,
                                                                       56128, 32948, 56143,
                                                                       14122, 14140, 33478,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 56968, 0, 3,
                                                                       56143, 32958, 56158,
                                                                       14140, 14158, 33508,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 57013, 0, 3,
                                                                       56158, 32968, 56173,
                                                                       14158, 14176, 33538,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 57058, 0, 3,
                                                                       56173, 32978, 56188,
                                                                       14176, 14194, 33568,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 57103, 0, 3,
                                                                       56188, 32988, 56203,
                                                                       14194, 14212, 33598,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 57148, 0, 3,
                                                                       56203, 32998, 56218,
                                                                       14212, 14230, 33628,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 57193, 0, 3,
                                                                       56218, 33008, 56233,
                                                                       14230, 14248, 33658,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 57238, 0, 3,
                                                                       56248, 33028, 56293,
                                                                       14284, 14320, 33688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 57328, 0, 3,
                                                                       56293, 33058, 56338,
                                                                       14320, 14356, 33748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 57418, 0, 3,
                                                                       56338, 33088, 56383,
                                                                       14356, 14392, 33808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 57508, 0, 3,
                                                                       56383, 33118, 56428,
                                                                       14392, 14428, 33868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 57598, 0, 3,
                                                                       56428, 33148, 56473,
                                                                       14428, 14464, 33928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 57688, 0, 3,
                                                                       56473, 33178, 56518,
                                                                       14464, 14500, 33988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 57778, 0, 3,
                                                                       56518, 33208, 56563,
                                                                       14500, 14536, 34048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 57868, 0, 3,
                                                                       56563, 33238, 56608,
                                                                       14536, 14572, 34108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 57958, 0, 3,
                                                                       56608, 33268, 56653,
                                                                       14572, 14608, 34168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 58048, 0, 3,
                                                                       56653, 33298, 56698,
                                                                       14608, 14644, 34228,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 58138, 0, 3,
                                                                       56743, 33358, 56788,
                                                                       14716, 14752, 34288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 58228, 0, 3,
                                                                       56788, 33388, 56833,
                                                                       14752, 14788, 34348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 58318, 0, 3,
                                                                       56833, 33418, 56878,
                                                                       14788, 14824, 34408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 58408, 0, 3,
                                                                       56878, 33448, 56923,
                                                                       14824, 14860, 34468,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 58498, 0, 3,
                                                                       56923, 33478, 56968,
                                                                       14860, 14896, 34528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 58588, 0, 3,
                                                                       56968, 33508, 57013,
                                                                       14896, 14932, 34588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 58678, 0, 3,
                                                                       57013, 33538, 57058,
                                                                       14932, 14968, 34648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 58768, 0, 3,
                                                                       57058, 33568, 57103,
                                                                       14968, 15004, 34708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 58858, 0, 3,
                                                                       57103, 33598, 57148,
                                                                       15004, 15040, 34768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 58948, 0, 3,
                                                                       57148, 33628, 57193,
                                                                       15040, 15076, 34828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 59038, 0, 3,
                                                                       57238, 33688, 57328,
                                                                       15148, 15208, 34888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 59188, 0, 3,
                                                                       57328, 33748, 57418,
                                                                       15208, 15268, 34988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 59338, 0, 3,
                                                                       57418, 33808, 57508,
                                                                       15268, 15328, 35088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 59488, 0, 3,
                                                                       57508, 33868, 57598,
                                                                       15328, 15388, 35188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 59638, 0, 3,
                                                                       57598, 33928, 57688,
                                                                       15388, 15448, 35288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 59788, 0, 3,
                                                                       57688, 33988, 57778,
                                                                       15448, 15508, 35388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 59938, 0, 3,
                                                                       57778, 34048, 57868,
                                                                       15508, 15568, 35488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 60088, 0, 3,
                                                                       57868, 34108, 57958,
                                                                       15568, 15628, 35588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 60238, 0, 3,
                                                                       57958, 34168, 58048,
                                                                       15628, 15688, 35688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 60388, 0, 3,
                                                                       58138, 34288, 58228,
                                                                       15808, 15868, 35788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 60538, 0, 3,
                                                                       58228, 34348, 58318,
                                                                       15868, 15928, 35888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 60688, 0, 3,
                                                                       58318, 34408, 58408,
                                                                       15928, 15988, 35988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 60838, 0, 3,
                                                                       58408, 34468, 58498,
                                                                       15988, 16048, 36088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 60988, 0, 3,
                                                                       58498, 34528, 58588,
                                                                       16048, 16108, 36188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61138, 0, 3,
                                                                       58588, 34588, 58678,
                                                                       16108, 16168, 36288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61288, 0, 3,
                                                                       58678, 34648, 58768,
                                                                       16168, 16228, 36388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61438, 0, 3,
                                                                       58768, 34708, 58858,
                                                                       16228, 16288, 36488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61588, 0, 3,
                                                                       58858, 34768, 58948,
                                                                       16288, 16348, 36588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 61738, 0, 3,
                                                                       59038, 34888, 59188,
                                                                       16468, 16558, 36688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 61963, 0, 3,
                                                                       59188, 34988, 59338,
                                                                       16558, 16648, 36838,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 62188, 0, 3,
                                                                       59338, 35088, 59488,
                                                                       16648, 16738, 36988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 62413, 0, 3,
                                                                       59488, 35188, 59638,
                                                                       16738, 16828, 37138,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 62638, 0, 3,
                                                                       59638, 35288, 59788,
                                                                       16828, 16918, 37288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 62863, 0, 3,
                                                                       59788, 35388, 59938,
                                                                       16918, 17008, 37438,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63088, 0, 3,
                                                                       59938, 35488, 60088,
                                                                       17008, 17098, 37588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63313, 0, 3,
                                                                       60088, 35588, 60238,
                                                                       17098, 17188, 37738,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63538, 0, 3,
                                                                       60388, 35788, 60538,
                                                                       17368, 17458, 37888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63763, 0, 3,
                                                                       60538, 35888, 60688,
                                                                       17458, 17548, 38038,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63988, 0, 3,
                                                                       60688, 35988, 60838,
                                                                       17548, 17638, 38188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 64213, 0, 3,
                                                                       60838, 36088, 60988,
                                                                       17638, 17728, 38338,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 64438, 0, 3,
                                                                       60988, 36188, 61138,
                                                                       17728, 17818, 38488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 64663, 0, 3,
                                                                       61138, 36288, 61288,
                                                                       17818, 17908, 38638,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 64888, 0, 3,
                                                                       61288, 36388, 61438,
                                                                       17908, 17998, 38788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 65113, 0, 3,
                                                                       61438, 36488, 61588,
                                                                       17998, 18088, 38938,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 65338, 0, 3,
                                                                       61738, 36688, 61963,
                                                                       18268, 18394, 39088,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 65653, 0, 3,
                                                                       61963, 36838, 62188,
                                                                       18394, 18520, 39298,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 65968, 0, 3,
                                                                       62188, 36988, 62413,
                                                                       18520, 18646, 39508,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 66283, 0, 3,
                                                                       62413, 37138, 62638,
                                                                       18646, 18772, 39718,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 66598, 0, 3,
                                                                       62638, 37288, 62863,
                                                                       18772, 18898, 39928,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 66913, 0, 3,
                                                                       62863, 37438, 63088,
                                                                       18898, 19024, 40138,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 67228, 0, 3,
                                                                       63088, 37588, 63313,
                                                                       19024, 19150, 40348,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 67543, 0, 3,
                                                                       63538, 37888, 63763,
                                                                       19402, 19528, 40558,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 67858, 0, 3,
                                                                       63763, 38038, 63988,
                                                                       19528, 19654, 40768,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 68173, 0, 3,
                                                                       63988, 38188, 64213,
                                                                       19654, 19780, 40978,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 68488, 0, 3,
                                                                       64213, 38338, 64438,
                                                                       19780, 19906, 41188,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 68803, 0, 3,
                                                                       64438, 38488, 64663,
                                                                       19906, 20032, 41398,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 69118, 0, 3,
                                                                       64663, 38638, 64888,
                                                                       20032, 20158, 41608,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 69433, 0, 3,
                                                                       64888, 38788, 65113,
                                                                       20158, 20284, 41818,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 69748, 0, 3,
                                                                       65338, 39088, 65653,
                                                                       20536, 20704, 42028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 70168, 0, 3,
                                                                       65653, 39298, 65968,
                                                                       20704, 20872, 42308,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 70588, 0, 3,
                                                                       65968, 39508, 66283,
                                                                       20872, 21040, 42588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 71008, 0, 3,
                                                                       66283, 39718, 66598,
                                                                       21040, 21208, 42868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 71428, 0, 3,
                                                                       66598, 39928, 66913,
                                                                       21208, 21376, 43148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 71848, 0, 3,
                                                                       66913, 40138, 67228,
                                                                       21376, 21544, 43428,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 72268, 0, 3,
                                                                       67543, 40558, 67858,
                                                                       21880, 22048, 43708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 72688, 0, 3,
                                                                       67858, 40768, 68173,
                                                                       22048, 22216, 43988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 73108, 0, 3,
                                                                       68173, 40978, 68488,
                                                                       22216, 22384, 44268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 73528, 0, 3,
                                                                       68488, 41188, 68803,
                                                                       22384, 22552, 44548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 73948, 0, 3,
                                                                       68803, 41398, 69118,
                                                                       22552, 22720, 44828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 74368, 0, 3,
                                                                       69118, 41608, 69433,
                                                                       22720, 22888, 45108,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 74788, 0, 3,
                                                                       69748, 42028, 70168,
                                                                       23224, 23440, 45388,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 75328, 0, 3,
                                                                       70168, 42308, 70588,
                                                                       23440, 23656, 45748,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 75868, 0, 3,
                                                                       70588, 42588, 71008,
                                                                       23656, 23872, 46108,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 76408, 0, 3,
                                                                       71008, 42868, 71428,
                                                                       23872, 24088, 46468,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 76948, 0, 3,
                                                                       71428, 43148, 71848,
                                                                       24088, 24304, 46828,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 77488, 0, 3,
                                                                       72268, 43708, 72688,
                                                                       24736, 24952, 47188,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 78028, 0, 3,
                                                                       72688, 43988, 73108,
                                                                       24952, 25168, 47548,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 78568, 0, 3,
                                                                       73108, 44268, 73528,
                                                                       25168, 25384, 47908,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 79108, 0, 3,
                                                                       73528, 44548, 73948,
                                                                       25384, 25600, 48268,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 79648, 0, 3,
                                                                       73948, 44828, 74368,
                                                                       25600, 25816, 48628,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 80188, 0, 3,
                                                                       74788, 45388, 75328,
                                                                       26248, 26518, 48988,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 80863, 0, 3,
                                                                       75328, 45748, 75868,
                                                                       26518, 26788, 49438,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 81538, 0, 3,
                                                                       75868, 46108, 76408,
                                                                       26788, 27058, 49888,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 82213, 0, 3,
                                                                       76408, 46468, 76948,
                                                                       27058, 27328, 50338,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 82888, 0, 3,
                                                                       77488, 47188, 78028,
                                                                       27868, 28138, 50788,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 83563, 0, 3,
                                                                       78028, 47548, 78568,
                                                                       28138, 28408, 51238,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 84238, 0, 3,
                                                                       78568, 47908, 79108,
                                                                       28408, 28678, 51688,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 84913, 0, 3,
                                                                       79108, 48268, 79648,
                                                                       28678, 28948, 52138,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 85588, 0, 3,
                                                                       80188, 48988, 80863,
                                                                       29488, 29818, 52588,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 86413, 0, 3,
                                                                       80863, 49438, 81538,
                                                                       29818, 30148, 53138,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 87238, 0, 3,
                                                                       81538, 49888, 82213,
                                                                       30148, 30478, 53688,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 88063, 0, 3,
                                                                       82888, 50788, 83563,
                                                                       31138, 31468, 54238,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 88888, 0, 3,
                                                                       83563, 51238, 84238,
                                                                       31468, 31798, 54788,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 89713, 0, 3,
                                                                       84238, 51688, 84913,
                                                                       31798, 32128, 55338,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90538, 3, 32788,
                                                                       32798, 55918, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90559, 3, 32798,
                                                                       32808, 55933, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90580, 3, 32808,
                                                                       32818, 55948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90601, 3, 32818,
                                                                       32828, 55963, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90622, 3, 32828,
                                                                       32838, 55978, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90643, 3, 32838,
                                                                       32848, 55993, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90664, 3, 32848,
                                                                       32858, 56008, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90685, 3, 32858,
                                                                       32868, 56023, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90706, 3, 32868,
                                                                       32878, 56038, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90727, 3, 32878,
                                                                       32888, 56053, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90748, 3, 32908,
                                                                       32918, 56098, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90769, 3, 32918,
                                                                       32928, 56113, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90790, 3, 32928,
                                                                       32938, 56128, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90811, 3, 32938,
                                                                       32948, 56143, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90832, 3, 32948,
                                                                       32958, 56158, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90853, 3, 32958,
                                                                       32968, 56173, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90874, 3, 32968,
                                                                       32978, 56188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90895, 3, 32978,
                                                                       32988, 56203, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90916, 3, 32988,
                                                                       32998, 56218, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90937, 3, 32998,
                                                                       33008, 56233, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 90958, 0, 3,
                                                                       90538, 55918, 90559,
                                                                       33028, 33058, 56338,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91021, 0, 3,
                                                                       90559, 55933, 90580,
                                                                       33058, 33088, 56383,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91084, 0, 3,
                                                                       90580, 55948, 90601,
                                                                       33088, 33118, 56428,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91147, 0, 3,
                                                                       90601, 55963, 90622,
                                                                       33118, 33148, 56473,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91210, 0, 3,
                                                                       90622, 55978, 90643,
                                                                       33148, 33178, 56518,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91273, 0, 3,
                                                                       90643, 55993, 90664,
                                                                       33178, 33208, 56563,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91336, 0, 3,
                                                                       90664, 56008, 90685,
                                                                       33208, 33238, 56608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91399, 0, 3,
                                                                       90685, 56023, 90706,
                                                                       33238, 33268, 56653,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91462, 0, 3,
                                                                       90706, 56038, 90727,
                                                                       33268, 33298, 56698,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91525, 0, 3,
                                                                       90748, 56098, 90769,
                                                                       33358, 33388, 56833,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91588, 0, 3,
                                                                       90769, 56113, 90790,
                                                                       33388, 33418, 56878,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91651, 0, 3,
                                                                       90790, 56128, 90811,
                                                                       33418, 33448, 56923,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91714, 0, 3,
                                                                       90811, 56143, 90832,
                                                                       33448, 33478, 56968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91777, 0, 3,
                                                                       90832, 56158, 90853,
                                                                       33478, 33508, 57013,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91840, 0, 3,
                                                                       90853, 56173, 90874,
                                                                       33508, 33538, 57058,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91903, 0, 3,
                                                                       90874, 56188, 90895,
                                                                       33538, 33568, 57103,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91966, 0, 3,
                                                                       90895, 56203, 90916,
                                                                       33568, 33598, 57148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 92029, 0, 3,
                                                                       90916, 56218, 90937,
                                                                       33598, 33628, 57193,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92092, 0, 3,
                                                                       90958, 56338, 91021,
                                                                       33688, 33748, 57418,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92218, 0, 3,
                                                                       91021, 56383, 91084,
                                                                       33748, 33808, 57508,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92344, 0, 3,
                                                                       91084, 56428, 91147,
                                                                       33808, 33868, 57598,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92470, 0, 3,
                                                                       91147, 56473, 91210,
                                                                       33868, 33928, 57688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92596, 0, 3,
                                                                       91210, 56518, 91273,
                                                                       33928, 33988, 57778,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92722, 0, 3,
                                                                       91273, 56563, 91336,
                                                                       33988, 34048, 57868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92848, 0, 3,
                                                                       91336, 56608, 91399,
                                                                       34048, 34108, 57958,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92974, 0, 3,
                                                                       91399, 56653, 91462,
                                                                       34108, 34168, 58048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 93100, 0, 3,
                                                                       91525, 56833, 91588,
                                                                       34288, 34348, 58318,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 93226, 0, 3,
                                                                       91588, 56878, 91651,
                                                                       34348, 34408, 58408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 93352, 0, 3,
                                                                       91651, 56923, 91714,
                                                                       34408, 34468, 58498,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 93478, 0, 3,
                                                                       91714, 56968, 91777,
                                                                       34468, 34528, 58588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 93604, 0, 3,
                                                                       91777, 57013, 91840,
                                                                       34528, 34588, 58678,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 93730, 0, 3,
                                                                       91840, 57058, 91903,
                                                                       34588, 34648, 58768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 93856, 0, 3,
                                                                       91903, 57103, 91966,
                                                                       34648, 34708, 58858,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 93982, 0, 3,
                                                                       91966, 57148, 92029,
                                                                       34708, 34768, 58948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 94108, 0, 3,
                                                                       92092, 57418, 92218,
                                                                       34888, 34988, 59338,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 94318, 0, 3,
                                                                       92218, 57508, 92344,
                                                                       34988, 35088, 59488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 94528, 0, 3,
                                                                       92344, 57598, 92470,
                                                                       35088, 35188, 59638,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 94738, 0, 3,
                                                                       92470, 57688, 92596,
                                                                       35188, 35288, 59788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 94948, 0, 3,
                                                                       92596, 57778, 92722,
                                                                       35288, 35388, 59938,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 95158, 0, 3,
                                                                       92722, 57868, 92848,
                                                                       35388, 35488, 60088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 95368, 0, 3,
                                                                       92848, 57958, 92974,
                                                                       35488, 35588, 60238,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 95578, 0, 3,
                                                                       93100, 58318, 93226,
                                                                       35788, 35888, 60688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 95788, 0, 3,
                                                                       93226, 58408, 93352,
                                                                       35888, 35988, 60838,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 95998, 0, 3,
                                                                       93352, 58498, 93478,
                                                                       35988, 36088, 60988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 96208, 0, 3,
                                                                       93478, 58588, 93604,
                                                                       36088, 36188, 61138,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 96418, 0, 3,
                                                                       93604, 58678, 93730,
                                                                       36188, 36288, 61288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 96628, 0, 3,
                                                                       93730, 58768, 93856,
                                                                       36288, 36388, 61438,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 96838, 0, 3,
                                                                       93856, 58858, 93982,
                                                                       36388, 36488, 61588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 97048, 0, 3,
                                                                       94108, 59338, 94318,
                                                                       36688, 36838, 62188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 97363, 0, 3,
                                                                       94318, 59488, 94528,
                                                                       36838, 36988, 62413,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 97678, 0, 3,
                                                                       94528, 59638, 94738,
                                                                       36988, 37138, 62638,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 97993, 0, 3,
                                                                       94738, 59788, 94948,
                                                                       37138, 37288, 62863,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 98308, 0, 3,
                                                                       94948, 59938, 95158,
                                                                       37288, 37438, 63088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 98623, 0, 3,
                                                                       95158, 60088, 95368,
                                                                       37438, 37588, 63313,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 98938, 0, 3,
                                                                       95578, 60688, 95788,
                                                                       37888, 38038, 63988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 99253, 0, 3,
                                                                       95788, 60838, 95998,
                                                                       38038, 38188, 64213,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 99568, 0, 3,
                                                                       95998, 60988, 96208,
                                                                       38188, 38338, 64438,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 99883, 0, 3,
                                                                       96208, 61138, 96418,
                                                                       38338, 38488, 64663,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 100198, 0, 3,
                                                                       96418, 61288, 96628,
                                                                       38488, 38638, 64888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 100513, 0, 3,
                                                                       96628, 61438, 96838,
                                                                       38638, 38788, 65113,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 100828, 0, 3,
                                                                       97048, 62188, 97363,
                                                                       39088, 39298, 65968,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 101269, 0, 3,
                                                                       97363, 62413, 97678,
                                                                       39298, 39508, 66283,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 101710, 0, 3,
                                                                       97678, 62638, 97993,
                                                                       39508, 39718, 66598,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 102151, 0, 3,
                                                                       97993, 62863, 98308,
                                                                       39718, 39928, 66913,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 102592, 0, 3,
                                                                       98308, 63088, 98623,
                                                                       39928, 40138, 67228,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 103033, 0, 3,
                                                                       98938, 63988, 99253,
                                                                       40558, 40768, 68173,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 103474, 0, 3,
                                                                       99253, 64213, 99568,
                                                                       40768, 40978, 68488,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 103915, 0, 3,
                                                                       99568, 64438, 99883,
                                                                       40978, 41188, 68803,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 104356, 0, 3,
                                                                       99883, 64663, 100198,
                                                                       41188, 41398, 69118,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 104797, 0, 3,
                                                                       100198, 64888, 100513,
                                                                       41398, 41608, 69433,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 105238, 0, 3,
                                                                       100828, 65968, 101269,
                                                                       42028, 42308, 70588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 105826, 0, 3,
                                                                       101269, 66283, 101710,
                                                                       42308, 42588, 71008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 106414, 0, 3,
                                                                       101710, 66598, 102151,
                                                                       42588, 42868, 71428,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 107002, 0, 3,
                                                                       102151, 66913, 102592,
                                                                       42868, 43148, 71848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 107590, 0, 3,
                                                                       103033, 68173, 103474,
                                                                       43708, 43988, 73108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 108178, 0, 3,
                                                                       103474, 68488, 103915,
                                                                       43988, 44268, 73528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 108766, 0, 3,
                                                                       103915, 68803, 104356,
                                                                       44268, 44548, 73948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 109354, 0, 3,
                                                                       104356, 69118, 104797,
                                                                       44548, 44828, 74368,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 109942, 0, 3,
                                                                       105238, 70588, 105826,
                                                                       45388, 45748, 75868,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 110698, 0, 3,
                                                                       105826, 71008, 106414,
                                                                       45748, 46108, 76408,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 111454, 0, 3,
                                                                       106414, 71428, 107002,
                                                                       46108, 46468, 76948,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 112210, 0, 3,
                                                                       107590, 73108, 108178,
                                                                       47188, 47548, 78568,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 112966, 0, 3,
                                                                       108178, 73528, 108766,
                                                                       47548, 47908, 79108,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 113722, 0, 3,
                                                                       108766, 73948, 109354,
                                                                       47908, 48268, 79648,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 114478, 0, 3,
                                                                       109942, 75868, 110698,
                                                                       48988, 49438, 81538,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 115423, 0, 3,
                                                                       110698, 76408, 111454,
                                                                       49438, 49888, 82213,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 116368, 0, 3,
                                                                       112210, 78568, 112966,
                                                                       50788, 51238, 84238,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 117313, 0, 3,
                                                                       112966, 79108, 113722,
                                                                       51238, 51688, 84913,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 118258, 0, 3,
                                                                       114478, 81538, 115423,
                                                                       52588, 53138, 87238,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 119413, 0, 3,
                                                                       116368, 84238, 117313,
                                                                       54238, 54788, 89713,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120568, 3, 55888,
                                                                       55903, 90538, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120596, 3, 55903,
                                                                       55918, 90559, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120624, 3, 55918,
                                                                       55933, 90580, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120652, 3, 55933,
                                                                       55948, 90601, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120680, 3, 55948,
                                                                       55963, 90622, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120708, 3, 55963,
                                                                       55978, 90643, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120736, 3, 55978,
                                                                       55993, 90664, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120764, 3, 55993,
                                                                       56008, 90685, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120792, 3, 56008,
                                                                       56023, 90706, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120820, 3, 56023,
                                                                       56038, 90727, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120848, 3, 56068,
                                                                       56083, 90748, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120876, 3, 56083,
                                                                       56098, 90769, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120904, 3, 56098,
                                                                       56113, 90790, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120932, 3, 56113,
                                                                       56128, 90811, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120960, 3, 56128,
                                                                       56143, 90832, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 120988, 3, 56143,
                                                                       56158, 90853, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 121016, 3, 56158,
                                                                       56173, 90874, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 121044, 3, 56173,
                                                                       56188, 90895, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 121072, 3, 56188,
                                                                       56203, 90916, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 121100, 3, 56203,
                                                                       56218, 90937, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 121128, 0, 3,
                                                                       120568, 90538, 120596,
                                                                       56248, 56293, 90958,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 121212, 0, 3,
                                                                       120596, 90559, 120624,
                                                                       56293, 56338, 91021,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 121296, 0, 3,
                                                                       120624, 90580, 120652,
                                                                       56338, 56383, 91084,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 121380, 0, 3,
                                                                       120652, 90601, 120680,
                                                                       56383, 56428, 91147,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 121464, 0, 3,
                                                                       120680, 90622, 120708,
                                                                       56428, 56473, 91210,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 121548, 0, 3,
                                                                       120708, 90643, 120736,
                                                                       56473, 56518, 91273,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 121632, 0, 3,
                                                                       120736, 90664, 120764,
                                                                       56518, 56563, 91336,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 121716, 0, 3,
                                                                       120764, 90685, 120792,
                                                                       56563, 56608, 91399,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 121800, 0, 3,
                                                                       120792, 90706, 120820,
                                                                       56608, 56653, 91462,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 121884, 0, 3,
                                                                       120848, 90748, 120876,
                                                                       56743, 56788, 91525,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 121968, 0, 3,
                                                                       120876, 90769, 120904,
                                                                       56788, 56833, 91588,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 122052, 0, 3,
                                                                       120904, 90790, 120932,
                                                                       56833, 56878, 91651,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 122136, 0, 3,
                                                                       120932, 90811, 120960,
                                                                       56878, 56923, 91714,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 122220, 0, 3,
                                                                       120960, 90832, 120988,
                                                                       56923, 56968, 91777,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 122304, 0, 3,
                                                                       120988, 90853, 121016,
                                                                       56968, 57013, 91840,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 122388, 0, 3,
                                                                       121016, 90874, 121044,
                                                                       57013, 57058, 91903,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 122472, 0, 3,
                                                                       121044, 90895, 121072,
                                                                       57058, 57103, 91966,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 122556, 0, 3,
                                                                       121072, 90916, 121100,
                                                                       57103, 57148, 92029,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 122640, 0, 3,
                                                                       121128, 90958, 121212,
                                                                       57238, 57328, 92092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 122808, 0, 3,
                                                                       121212, 91021, 121296,
                                                                       57328, 57418, 92218,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 122976, 0, 3,
                                                                       121296, 91084, 121380,
                                                                       57418, 57508, 92344,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 123144, 0, 3,
                                                                       121380, 91147, 121464,
                                                                       57508, 57598, 92470,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 123312, 0, 3,
                                                                       121464, 91210, 121548,
                                                                       57598, 57688, 92596,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 123480, 0, 3,
                                                                       121548, 91273, 121632,
                                                                       57688, 57778, 92722,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 123648, 0, 3,
                                                                       121632, 91336, 121716,
                                                                       57778, 57868, 92848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 123816, 0, 3,
                                                                       121716, 91399, 121800,
                                                                       57868, 57958, 92974,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 123984, 0, 3,
                                                                       121884, 91525, 121968,
                                                                       58138, 58228, 93100,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 124152, 0, 3,
                                                                       121968, 91588, 122052,
                                                                       58228, 58318, 93226,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 124320, 0, 3,
                                                                       122052, 91651, 122136,
                                                                       58318, 58408, 93352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 124488, 0, 3,
                                                                       122136, 91714, 122220,
                                                                       58408, 58498, 93478,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 124656, 0, 3,
                                                                       122220, 91777, 122304,
                                                                       58498, 58588, 93604,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 124824, 0, 3,
                                                                       122304, 91840, 122388,
                                                                       58588, 58678, 93730,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 124992, 0, 3,
                                                                       122388, 91903, 122472,
                                                                       58678, 58768, 93856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 125160, 0, 3,
                                                                       122472, 91966, 122556,
                                                                       58768, 58858, 93982,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 125328, 0, 3,
                                                                       122640, 92092, 122808,
                                                                       59038, 59188, 94108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 125608, 0, 3,
                                                                       122808, 92218, 122976,
                                                                       59188, 59338, 94318,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 125888, 0, 3,
                                                                       122976, 92344, 123144,
                                                                       59338, 59488, 94528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 126168, 0, 3,
                                                                       123144, 92470, 123312,
                                                                       59488, 59638, 94738,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 126448, 0, 3,
                                                                       123312, 92596, 123480,
                                                                       59638, 59788, 94948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 126728, 0, 3,
                                                                       123480, 92722, 123648,
                                                                       59788, 59938, 95158,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 127008, 0, 3,
                                                                       123648, 92848, 123816,
                                                                       59938, 60088, 95368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 127288, 0, 3,
                                                                       123984, 93100, 124152,
                                                                       60388, 60538, 95578,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 127568, 0, 3,
                                                                       124152, 93226, 124320,
                                                                       60538, 60688, 95788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 127848, 0, 3,
                                                                       124320, 93352, 124488,
                                                                       60688, 60838, 95998,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 128128, 0, 3,
                                                                       124488, 93478, 124656,
                                                                       60838, 60988, 96208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 128408, 0, 3,
                                                                       124656, 93604, 124824,
                                                                       60988, 61138, 96418,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 128688, 0, 3,
                                                                       124824, 93730, 124992,
                                                                       61138, 61288, 96628,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 128968, 0, 3,
                                                                       124992, 93856, 125160,
                                                                       61288, 61438, 96838,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 129248, 0, 3,
                                                                       125328, 94108, 125608,
                                                                       61738, 61963, 97048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 129668, 0, 3,
                                                                       125608, 94318, 125888,
                                                                       61963, 62188, 97363,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 130088, 0, 3,
                                                                       125888, 94528, 126168,
                                                                       62188, 62413, 97678,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 130508, 0, 3,
                                                                       126168, 94738, 126448,
                                                                       62413, 62638, 97993,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 130928, 0, 3,
                                                                       126448, 94948, 126728,
                                                                       62638, 62863, 98308,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 131348, 0, 3,
                                                                       126728, 95158, 127008,
                                                                       62863, 63088, 98623,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 131768, 0, 3,
                                                                       127288, 95578, 127568,
                                                                       63538, 63763, 98938,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 132188, 0, 3,
                                                                       127568, 95788, 127848,
                                                                       63763, 63988, 99253,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 132608, 0, 3,
                                                                       127848, 95998, 128128,
                                                                       63988, 64213, 99568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 133028, 0, 3,
                                                                       128128, 96208, 128408,
                                                                       64213, 64438, 99883,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 133448, 0, 3,
                                                                       128408, 96418, 128688,
                                                                       64438, 64663, 100198,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 133868, 0, 3,
                                                                       128688, 96628, 128968,
                                                                       64663, 64888, 100513,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 134288, 0, 3,
                                                                       129248, 97048, 129668,
                                                                       65338, 65653, 100828,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 134876, 0, 3,
                                                                       129668, 97363, 130088,
                                                                       65653, 65968, 101269,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 135464, 0, 3,
                                                                       130088, 97678, 130508,
                                                                       65968, 66283, 101710,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 136052, 0, 3,
                                                                       130508, 97993, 130928,
                                                                       66283, 66598, 102151,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 136640, 0, 3,
                                                                       130928, 98308, 131348,
                                                                       66598, 66913, 102592,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 137228, 0, 3,
                                                                       131768, 98938, 132188,
                                                                       67543, 67858, 103033,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 137816, 0, 3,
                                                                       132188, 99253, 132608,
                                                                       67858, 68173, 103474,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 138404, 0, 3,
                                                                       132608, 99568, 133028,
                                                                       68173, 68488, 103915,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 138992, 0, 3,
                                                                       133028, 99883, 133448,
                                                                       68488, 68803, 104356,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 139580, 0, 3,
                                                                       133448, 100198, 133868,
                                                                       68803, 69118, 104797,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 140168, 0, 3,
                                                                       134288, 100828, 134876,
                                                                       69748, 70168, 105238,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 140952, 0, 3,
                                                                       134876, 101269, 135464,
                                                                       70168, 70588, 105826,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 141736, 0, 3,
                                                                       135464, 101710, 136052,
                                                                       70588, 71008, 106414,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 142520, 0, 3,
                                                                       136052, 102151, 136640,
                                                                       71008, 71428, 107002,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 143304, 0, 3,
                                                                       137228, 103033, 137816,
                                                                       72268, 72688, 107590,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 144088, 0, 3,
                                                                       137816, 103474, 138404,
                                                                       72688, 73108, 108178,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 144872, 0, 3,
                                                                       138404, 103915, 138992,
                                                                       73108, 73528, 108766,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 145656, 0, 3,
                                                                       138992, 104356, 139580,
                                                                       73528, 73948, 109354,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 146440, 0, 3,
                                                                       140168, 105238, 140952,
                                                                       74788, 75328, 109942,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 147448, 0, 3,
                                                                       140952, 105826, 141736,
                                                                       75328, 75868, 110698,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 148456, 0, 3,
                                                                       141736, 106414, 142520,
                                                                       75868, 76408, 111454,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 149464, 0, 3,
                                                                       143304, 107590, 144088,
                                                                       77488, 78028, 112210,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 150472, 0, 3,
                                                                       144088, 108178, 144872,
                                                                       78028, 78568, 112966,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 151480, 0, 3,
                                                                       144872, 108766, 145656,
                                                                       78568, 79108, 113722,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 152488, 0, 3,
                                                                       146440, 109942, 147448,
                                                                       80188, 80863, 114478,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 153748, 0, 3,
                                                                       147448, 110698, 148456,
                                                                       80863, 81538, 115423,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 155008, 0, 3,
                                                                       149464, 112210, 150472,
                                                                       82888, 83563, 116368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 156268, 0, 3,
                                                                       150472, 112966, 151480,
                                                                       83563, 84238, 117313,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 157528, 0, 3,
                                                                       152488, 114478, 153748,
                                                                       85588, 86413, 118258,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 159068, 0, 3,
                                                                       155008, 116368, 156268,
                                                                       88063, 88888, 119413,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_g_x(buffer, 160608, 125328, 134288, 1, 28, ncols, beta);

                    simdgeo::geom_g_y(buffer, 161028, 125328, 134288, 1, 28, ncols, beta);

                    simdgeo::geom_g_z(buffer, 161448, 125328, 134288, 1, 28, ncols, beta);

                    simdgeo::geom_g_x(buffer, 161868, 127288, 137228, 1, 28, ncols, beta);

                    simdgeo::geom_g_y(buffer, 162288, 127288, 137228, 1, 28, ncols, beta);

                    simdgeo::geom_g_z(buffer, 162708, 127288, 137228, 1, 28, ncols, beta);

                    simdgeo::geom_h_x(buffer, 163128, 129248, 140168, 1, 28, ncols, beta);

                    simdgeo::geom_h_y(buffer, 163716, 129248, 140168, 1, 28, ncols, beta);

                    simdgeo::geom_h_z(buffer, 164304, 129248, 140168, 1, 28, ncols, beta);

                    simdgeo::geom_h_x(buffer, 164892, 131768, 143304, 1, 28, ncols, beta);

                    simdgeo::geom_h_y(buffer, 165480, 131768, 143304, 1, 28, ncols, beta);

                    simdgeo::geom_h_z(buffer, 166068, 131768, 143304, 1, 28, ncols, beta);

                    simdgeo::geom_i_x(buffer, 166656, 134288, 146440, 1, 28, ncols, beta);

                    simdgeo::geom_i_y(buffer, 167440, 134288, 146440, 1, 28, ncols, beta);

                    simdgeo::geom_i_z(buffer, 168224, 134288, 146440, 1, 28, ncols, beta);

                    simdgeo::geom_i_x(buffer, 169008, 137228, 149464, 1, 28, ncols, beta);

                    simdgeo::geom_i_y(buffer, 169792, 137228, 149464, 1, 28, ncols, beta);

                    simdgeo::geom_i_z(buffer, 170576, 137228, 149464, 1, 28, ncols, beta);

                    simdgeo::geom_k_x(buffer, 171360, 140168, 152488, 1, 28, ncols, beta);

                    simdgeo::geom_k_y(buffer, 172368, 140168, 152488, 1, 28, ncols, beta);

                    simdgeo::geom_k_z(buffer, 173376, 140168, 152488, 1, 28, ncols, beta);

                    simdgeo::geom_k_x(buffer, 174384, 143304, 155008, 1, 28, ncols, beta);

                    simdgeo::geom_k_y(buffer, 175392, 143304, 155008, 1, 28, ncols, beta);

                    simdgeo::geom_k_z(buffer, 176400, 143304, 155008, 1, 28, ncols, beta);

                    simdgeo::geom_l_x(buffer, 177408, 146440, 157528, 1, 28, ncols, beta);

                    simdgeo::geom_l_y(buffer, 178668, 146440, 157528, 1, 28, ncols, beta);

                    simdgeo::geom_l_z(buffer, 179928, 146440, 157528, 1, 28, ncols, beta);

                    simdgeo::geom_l_x(buffer, 181188, 149464, 159068, 1, 28, ncols, beta);

                    simdgeo::geom_l_y(buffer, 182448, 149464, 159068, 1, 28, ncols, beta);

                    simdgeo::geom_l_z(buffer, 183708, 149464, 159068, 1, 28, ncols, beta);

                    simdfunc::contract_primitives(buffer, 184968, 160608, 420, ncols);

                    simdfunc::contract_primitives(buffer, 185583, 161028, 420, ncols);

                    simdfunc::contract_primitives(buffer, 186198, 161448, 420, ncols);

                    simdfunc::contract_primitives(buffer, 186813, 129248, 420, ncols);

                    simdfunc::contract_primitives(buffer, 187428, 161868, 420, ncols);

                    simdfunc::contract_primitives(buffer, 188043, 162288, 420, ncols);

                    simdfunc::contract_primitives(buffer, 188658, 162708, 420, ncols);

                    simdfunc::contract_primitives(buffer, 189273, 131768, 420, ncols);

                    simdfunc::contract_primitives(buffer, 189888, 163128, 588, ncols);

                    simdfunc::contract_primitives(buffer, 190749, 163716, 588, ncols);

                    simdfunc::contract_primitives(buffer, 191610, 164304, 588, ncols);

                    simdfunc::contract_primitives(buffer, 192471, 134288, 588, ncols);

                    simdfunc::contract_primitives(buffer, 193332, 164892, 588, ncols);

                    simdfunc::contract_primitives(buffer, 194193, 165480, 588, ncols);

                    simdfunc::contract_primitives(buffer, 195054, 166068, 588, ncols);

                    simdfunc::contract_primitives(buffer, 195915, 137228, 588, ncols);

                    simdfunc::contract_primitives(buffer, 196776, 166656, 784, ncols);

                    simdfunc::contract_primitives(buffer, 197924, 167440, 784, ncols);

                    simdfunc::contract_primitives(buffer, 199072, 168224, 784, ncols);

                    simdfunc::contract_primitives(buffer, 200220, 140168, 784, ncols);

                    simdfunc::contract_primitives(buffer, 201368, 169008, 784, ncols);

                    simdfunc::contract_primitives(buffer, 202516, 169792, 784, ncols);

                    simdfunc::contract_primitives(buffer, 203664, 170576, 784, ncols);

                    simdfunc::contract_primitives(buffer, 204812, 143304, 784, ncols);

                    simdfunc::contract_primitives(buffer, 205960, 171360, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 207436, 172368, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 208912, 173376, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 210388, 146440, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 211864, 174384, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 213340, 175392, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 214816, 176400, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 216292, 149464, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 217768, 177408, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 219613, 178668, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 221458, 179928, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 223303, 181188, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 225148, 182448, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 226993, 183708, 1260, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 185388, 184968, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 186003, 185583, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 186618, 186198, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 187233, 186813, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 187848, 187428, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 188463, 188043, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 189078, 188658, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 189693, 189273, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 190476, 189888, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 191337, 190749, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 192198, 191610, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 193059, 192471, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 193920, 193332, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 194781, 194193, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 195642, 195054, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 196503, 195915, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 197560, 196776, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 198708, 197924, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 199856, 199072, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 201004, 200220, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 202152, 201368, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 203300, 202516, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 204448, 203664, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 205596, 204812, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 206968, 205960, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 208444, 207436, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 209920, 208912, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 211396, 210388, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 212872, 211864, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 214348, 213340, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 215824, 214816, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 217300, 216292, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 219028, 217768, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 220873, 219613, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 222718, 221458, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 224563, 223303, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 226408, 225148, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 228253, 226993, 45, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 228838, 185388, 187233, 190476,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 229423, 186003, 187233, 191337,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 230008, 186618, 187233, 192198,
                                          13, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 230593, 187233, 193059, 13, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 231178, 187848, 189693, 193920,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 231763, 188463, 189693, 194781,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 232348, 189078, 189693, 195642,
                                          13, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 232933, 189693, 196503, 13, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 233518, 190476, 193059, 197560,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 234337, 191337, 193059, 198708,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 235156, 192198, 193059, 199856,
                                          13, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 235975, 193059, 201004, 13, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 236794, 193920, 196503, 202152,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 237613, 194781, 196503, 203300,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 238432, 195642, 196503, 204448,
                                          13, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 239251, 196503, 205596, 13, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 240070, 197560, 201004, 206968,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 241162, 198708, 201004, 208444,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 242254, 199856, 201004, 209920,
                                          13, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 243346, 201004, 211396, 13, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 244438, 202152, 205596, 212872,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 245530, 203300, 205596, 214348,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 246622, 204448, 205596, 215824,
                                          13, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 247714, 205596, 217300, 13, nmax);

        simdtrf::compute_hrr_geom_010x_pk(buffer, coordinates, 248806, 206968, 211396, 219028,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_pk(buffer, coordinates, 250210, 208444, 211396, 220873,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_pk(buffer, coordinates, 251614, 209920, 211396, 222718,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010x_pk(buffer, coordinates, 253018, 212872, 217300, 224563,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_pk(buffer, coordinates, 254422, 214348, 217300, 226408,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_pk(buffer, coordinates, 255826, 215824, 217300, 228253,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 257230, 228838, 230593, 233518,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 258400, 229423, 230593, 234337,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 259570, 230008, 230593, 235156,
                                          13, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 260740, 230593, 235975, 13, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 261910, 231178, 232933, 236794,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 263080, 231763, 232933, 237613,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 264250, 232348, 232933, 238432,
                                          13, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 265420, 232933, 239251, 13, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 266590, 233518, 235975, 240070,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 268228, 234337, 235975, 241162,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 269866, 235156, 235975, 242254,
                                          13, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 271504, 235975, 243346, 13, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 273142, 236794, 239251, 244438,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 274780, 237613, 239251, 245530,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 276418, 238432, 239251, 246622,
                                          13, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 278056, 239251, 247714, 13, nmax);

        simdtrf::compute_hrr_geom_010x_di(buffer, coordinates, 279694, 240070, 243346, 248806,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_di(buffer, coordinates, 281878, 241162, 243346, 250210,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_di(buffer, coordinates, 284062, 242254, 243346, 251614,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010x_di(buffer, coordinates, 286246, 244438, 247714, 253018,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_di(buffer, coordinates, 288430, 245530, 247714, 254422,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_di(buffer, coordinates, 290614, 246622, 247714, 255826,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 292798, 257230, 260740, 266590,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 294748, 258400, 260740, 268228,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 296698, 259570, 260740, 269866,
                                          13, nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 298648, 260740, 271504, 13, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 300598, 261910, 265420, 273142,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 302548, 263080, 265420, 274780,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 304498, 264250, 265420, 276418,
                                          13, nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 306448, 265420, 278056, 13, nmax);

        simdtrf::compute_hrr_geom_010x_fh(buffer, coordinates, 308398, 266590, 271504, 279694,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_fh(buffer, coordinates, 311128, 268228, 271504, 281878,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_fh(buffer, coordinates, 313858, 269866, 271504, 284062,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010x_fh(buffer, coordinates, 316588, 273142, 278056, 286246,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_fh(buffer, coordinates, 319318, 274780, 278056, 288430,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_fh(buffer, coordinates, 322048, 276418, 278056, 290614,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010x_gg(buffer, coordinates, 324778, 292798, 298648, 308398,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_gg(buffer, coordinates, 327703, 294748, 298648, 311128,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_gg(buffer, coordinates, 330628, 296698, 298648, 313858,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010x_gg(buffer, coordinates, 333553, 300598, 306448, 316588,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010y_gg(buffer, coordinates, 336478, 302548, 306448, 319318,
                                          13, nmax);

        simdtrf::compute_hrr_geom_010z_gg(buffer, coordinates, 339403, 304498, 306448, 322048,
                                          13, nmax);

        simdtrf::transform_g_inner(buffer, 342328, 333553, 15, 13, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 342328, 117, nmax);

        simdtrf::transform_g_inner(buffer, 342328, 336478, 15, 13, nmax);

        simdtrf::transform_g_outer(values + 1053 * nvalues + n * npairs, nvalues, buffer, 342328,
                                   117, nmax);

        simdtrf::transform_g_inner(buffer, 342328, 339403, 15, 13, nmax);

        simdtrf::transform_g_outer(values + 2106 * nvalues + n * npairs, nvalues, buffer, 342328,
                                   117, nmax);

        simdtrf::transform_g_inner(buffer, 342328, 324778, 15, 13, nmax);

        simdtrf::transform_g_outer(values + 3159 * nvalues + n * npairs, nvalues, buffer, 342328,
                                   117, nmax);

        simdtrf::transform_g_inner(buffer, 342328, 327703, 15, 13, nmax);

        simdtrf::transform_g_outer(values + 4212 * nvalues + n * npairs, nvalues, buffer, 342328,
                                   117, nmax);

        simdtrf::transform_g_inner(buffer, 342328, 330628, 15, 13, nmax);

        simdtrf::transform_g_outer(values + 5265 * nvalues + n * npairs, nvalues, buffer, 342328,
                                   117, nmax);
    }

    for (size_t m = 0; m < 6318; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
