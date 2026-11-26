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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecFGH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferFD.hpp"
#include "SimdTransferFF.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGD.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferGeom100XFD.hpp"
#include "SimdTransferGeom100XFF.hpp"
#include "SimdTransferGeom100XFG.hpp"
#include "SimdTransferGeom100XFP.hpp"
#include "SimdTransferGeom100XGD.hpp"
#include "SimdTransferGeom100XGF.hpp"
#include "SimdTransferGeom100XGP.hpp"
#include "SimdTransferGeom100XHD.hpp"
#include "SimdTransferGeom100XHP.hpp"
#include "SimdTransferGeom100XIP.hpp"
#include "SimdTransferGeom100YFD.hpp"
#include "SimdTransferGeom100YFF.hpp"
#include "SimdTransferGeom100YFG.hpp"
#include "SimdTransferGeom100YFP.hpp"
#include "SimdTransferGeom100YGD.hpp"
#include "SimdTransferGeom100YGF.hpp"
#include "SimdTransferGeom100YGP.hpp"
#include "SimdTransferGeom100YHD.hpp"
#include "SimdTransferGeom100YHP.hpp"
#include "SimdTransferGeom100YIP.hpp"
#include "SimdTransferGeom100ZFD.hpp"
#include "SimdTransferGeom100ZFF.hpp"
#include "SimdTransferGeom100ZFG.hpp"
#include "SimdTransferGeom100ZFP.hpp"
#include "SimdTransferGeom100ZGD.hpp"
#include "SimdTransferGeom100ZGF.hpp"
#include "SimdTransferGeom100ZGP.hpp"
#include "SimdTransferGeom100ZHD.hpp"
#include "SimdTransferGeom100ZHP.hpp"
#include "SimdTransferGeom100ZIP.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_fgh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_fgh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 179798, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 179798, 85148, 25460, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto alpha = a_exps[i];

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fa = -b_exps[j] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pa(buffer, coordinates, 0, nmax, fa);

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

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 76, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 79, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 85, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 88, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 91, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 94, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 97, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 100, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 106, 0, 3, 7, 8,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 8, 9,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 9, 10,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 124, 0, 3, 10, 11,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 130, 0, 3, 11, 12,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 136, 0, 3, 12, 13,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 13, 14,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 14, 15,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 154, 0, 3, 15, 16,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 160, 0, 3, 16, 17,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 166, 0, 3, 17, 18,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 21, 22,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 178, 0, 3, 22, 23,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 184, 0, 3, 23, 24,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 190, 0, 3, 24, 25,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 196, 0, 3, 25, 26,
                                                                       82, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 26, 27,
                                                                       85, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 208, 0, 3, 27, 28,
                                                                       88, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 214, 0, 3, 28, 29,
                                                                       91, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 220, 0, 3, 29, 30,
                                                                       94, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 226, 0, 3, 30, 31,
                                                                       97, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 31, 32,
                                                                       100, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 34, 37,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 37, 40,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 40, 43,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 43, 46,
                                                                       124, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 46, 49,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 49, 52,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 52, 55,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 55, 58,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 58, 61,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 61, 64,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 70, 73,
                                                                       172, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 73, 76,
                                                                       178, 184, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 76, 79,
                                                                       184, 190, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 79, 82,
                                                                       190, 196, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 82, 85,
                                                                       196, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 85, 88,
                                                                       202, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 88, 91,
                                                                       208, 214, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 91, 94,
                                                                       214, 220, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 94, 97,
                                                                       220, 226, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 97,
                                                                       100, 226, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 106,
                                                                       112, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 453, 0, 3, 112,
                                                                       118, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 118,
                                                                       124, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 483, 0, 3, 124,
                                                                       130, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 130,
                                                                       136, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 513, 0, 3, 136,
                                                                       142, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 528, 0, 3, 142,
                                                                       148, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 543, 0, 3, 148,
                                                                       154, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 558, 0, 3, 154,
                                                                       160, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 573, 0, 3, 172,
                                                                       178, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 588, 0, 3, 178,
                                                                       184, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 603, 0, 3, 184,
                                                                       190, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 618, 0, 3, 190,
                                                                       196, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 633, 0, 3, 196,
                                                                       202, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 648, 0, 3, 202,
                                                                       208, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 663, 0, 3, 208,
                                                                       214, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 678, 0, 3, 214,
                                                                       220, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 693, 0, 3, 220,
                                                                       226, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 708, 0, 3, 238,
                                                                       248, 438, 453, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 729, 0, 3, 248,
                                                                       258, 453, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 750, 0, 3, 258,
                                                                       268, 468, 483, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 771, 0, 3, 268,
                                                                       278, 483, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 792, 0, 3, 278,
                                                                       288, 498, 513, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 813, 0, 3, 288,
                                                                       298, 513, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 834, 0, 3, 298,
                                                                       308, 528, 543, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 855, 0, 3, 308,
                                                                       318, 543, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 876, 0, 3, 338,
                                                                       348, 573, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 897, 0, 3, 348,
                                                                       358, 588, 603, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 918, 0, 3, 358,
                                                                       368, 603, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 939, 0, 3, 368,
                                                                       378, 618, 633, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 960, 0, 3, 378,
                                                                       388, 633, 648, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 981, 0, 3, 388,
                                                                       398, 648, 663, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 398,
                                                                       408, 663, 678, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1023, 0, 3, 408,
                                                                       418, 678, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 438,
                                                                       453, 708, 729, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 453,
                                                                       468, 729, 750, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 468,
                                                                       483, 750, 771, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 483,
                                                                       498, 771, 792, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 498,
                                                                       513, 792, 813, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 513,
                                                                       528, 813, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 528,
                                                                       543, 834, 855, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 573,
                                                                       588, 876, 897, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 588,
                                                                       603, 897, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 603,
                                                                       618, 918, 939, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 618,
                                                                       633, 939, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 633,
                                                                       648, 960, 981, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 648,
                                                                       663, 981, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 663,
                                                                       678, 1002, 1023, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 708,
                                                                       729, 1044, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1472, 0, 3, 729,
                                                                       750, 1072, 1100, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1508, 0, 3, 750,
                                                                       771, 1100, 1128, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1544, 0, 3, 771,
                                                                       792, 1128, 1156, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1580, 0, 3, 792,
                                                                       813, 1156, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1616, 0, 3, 813,
                                                                       834, 1184, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1652, 0, 3, 876,
                                                                       897, 1240, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 897,
                                                                       918, 1268, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1724, 0, 3, 918,
                                                                       939, 1296, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1760, 0, 3, 939,
                                                                       960, 1324, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1796, 0, 3, 960,
                                                                       981, 1352, 1380, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1832, 0, 3, 981,
                                                                       1002, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1868, 0, 3, 1044,
                                                                       1072, 1436, 1472, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1913, 0, 3, 1072,
                                                                       1100, 1472, 1508, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1958, 0, 3, 1100,
                                                                       1128, 1508, 1544, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2003, 0, 3, 1128,
                                                                       1156, 1544, 1580, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2048, 0, 3, 1156,
                                                                       1184, 1580, 1616, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2093, 0, 3, 1240,
                                                                       1268, 1652, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2138, 0, 3, 1268,
                                                                       1296, 1688, 1724, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2183, 0, 3, 1296,
                                                                       1324, 1724, 1760, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2228, 0, 3, 1324,
                                                                       1352, 1760, 1796, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2273, 0, 3, 1352,
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

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2396, 3, 7, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2405, 3, 8, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2414, 3, 9, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2423, 3, 10, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2432, 3, 11, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2441, 3, 12, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2450, 3, 13, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2459, 3, 14, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2468, 3, 15, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2477, 3, 16, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2486, 3, 17, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2495, 3, 18, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2504, 3, 21, 70,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2513, 3, 22, 73,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2522, 3, 23, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2531, 3, 24, 79,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2540, 3, 25, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2549, 3, 26, 85,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2558, 3, 27, 88,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2567, 3, 28, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2576, 3, 29, 94,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2585, 3, 30, 97,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2594, 3, 31, 100,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2603, 3, 32, 103,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2612, 3, 34, 106,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2630, 3, 37, 112,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2648, 3, 40, 118,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2666, 3, 43, 124,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2684, 3, 46, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2702, 3, 49, 136,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2720, 3, 52, 142,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2738, 3, 55, 148,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2756, 3, 58, 154,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2774, 3, 61, 160,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2792, 3, 64, 166,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2810, 3, 70, 172,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2828, 3, 73, 178,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2846, 3, 76, 184,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2864, 3, 79, 190,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2882, 3, 82, 196,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2900, 3, 85, 202,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2918, 3, 88, 208,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2936, 3, 91, 214,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2954, 3, 94, 220,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2972, 3, 97, 226,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2990, 3, 100, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3008, 3, 106, 238,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3038, 3, 112, 248,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3068, 3, 118, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3098, 3, 124, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3128, 3, 130, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3158, 3, 136, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3188, 3, 142, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3218, 3, 148, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3248, 3, 154, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3278, 3, 160, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3308, 3, 172, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3338, 3, 178, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3368, 3, 184, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3398, 3, 190, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3428, 3, 196, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3458, 3, 202, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3488, 3, 208, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3518, 3, 214, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3548, 3, 220, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3578, 3, 226, 428,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3608, 3, 238, 438,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3653, 3, 248, 453,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3698, 3, 258, 468,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3743, 3, 268, 483,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3788, 3, 278, 498,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3833, 3, 288, 513,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3878, 3, 298, 528,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3923, 3, 308, 543,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3968, 3, 318, 558,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4013, 3, 338, 573,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4058, 3, 348, 588,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4103, 3, 358, 603,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4148, 3, 368, 618,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4193, 3, 378, 633,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4238, 3, 388, 648,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4283, 3, 398, 663,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4328, 3, 408, 678,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4373, 3, 418, 693,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4418, 3, 438, 708,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4481, 3, 453, 729,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4544, 3, 468, 750,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4607, 3, 483, 771,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4670, 3, 498, 792,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4733, 3, 513, 813,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4796, 3, 528, 834,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4859, 3, 543, 855,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4922, 3, 573, 876,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4985, 3, 588, 897,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5048, 3, 603, 918,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5111, 3, 618, 939,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5174, 3, 633, 960,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5237, 3, 648, 981,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5300, 3, 663,
                                                                       1002, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5363, 3, 678,
                                                                       1023, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5426, 3, 708,
                                                                       1044, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5510, 3, 729,
                                                                       1072, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5594, 3, 750,
                                                                       1100, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5678, 3, 771,
                                                                       1128, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5762, 3, 792,
                                                                       1156, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5846, 3, 813,
                                                                       1184, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5930, 3, 834,
                                                                       1212, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6014, 3, 876,
                                                                       1240, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6098, 3, 897,
                                                                       1268, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6182, 3, 918,
                                                                       1296, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6266, 3, 939,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6350, 3, 960,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6434, 3, 981,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6518, 3, 1002,
                                                                       1408, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6602, 3, 1044,
                                                                       1436, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6710, 3, 1072,
                                                                       1472, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6818, 3, 1100,
                                                                       1508, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6926, 3, 1128,
                                                                       1544, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7034, 3, 1156,
                                                                       1580, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7142, 3, 1184,
                                                                       1616, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7250, 3, 1240,
                                                                       1652, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7358, 3, 1268,
                                                                       1688, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7466, 3, 1296,
                                                                       1724, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7574, 3, 1324,
                                                                       1760, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7682, 3, 1352,
                                                                       1796, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7790, 3, 1380,
                                                                       1832, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7898, 3, 1436,
                                                                       1868, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8033, 3, 1472,
                                                                       1913, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8168, 3, 1508,
                                                                       1958, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8303, 3, 1544,
                                                                       2003, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8438, 3, 1580,
                                                                       2048, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8573, 3, 1652,
                                                                       2093, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8708, 3, 1688,
                                                                       2138, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8843, 3, 1724,
                                                                       2183, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8978, 3, 1760,
                                                                       2228, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9113, 3, 1796,
                                                                       2273, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9248, 3, 7, 8,
                                                                       2324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9254, 3, 8, 9,
                                                                       2327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9260, 3, 9, 10,
                                                                       2330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9266, 3, 10, 11,
                                                                       2333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9272, 3, 11, 12,
                                                                       2336, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9278, 3, 12, 13,
                                                                       2339, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9284, 3, 13, 14,
                                                                       2342, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9290, 3, 14, 15,
                                                                       2345, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9296, 3, 15, 16,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9302, 3, 16, 17,
                                                                       2351, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9308, 3, 17, 18,
                                                                       2354, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9314, 3, 21, 22,
                                                                       2363, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9320, 3, 22, 23,
                                                                       2366, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9326, 3, 23, 24,
                                                                       2369, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9332, 3, 24, 25,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9338, 3, 25, 26,
                                                                       2375, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9344, 3, 26, 27,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9350, 3, 27, 28,
                                                                       2381, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9356, 3, 28, 29,
                                                                       2384, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9362, 3, 29, 30,
                                                                       2387, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9368, 3, 30, 31,
                                                                       2390, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9374, 3, 31, 32,
                                                                       2393, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9380, 0, 3, 9248,
                                                                       2324, 9254, 2414, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9398, 0, 3, 9254,
                                                                       2327, 9260, 2423, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9416, 0, 3, 9260,
                                                                       2330, 9266, 2432, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9434, 0, 3, 9266,
                                                                       2333, 9272, 2441, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9452, 0, 3, 9272,
                                                                       2336, 9278, 2450, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9470, 0, 3, 9278,
                                                                       2339, 9284, 2459, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9488, 0, 3, 9284,
                                                                       2342, 9290, 2468, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9506, 0, 3, 9290,
                                                                       2345, 9296, 2477, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9524, 0, 3, 9296,
                                                                       2348, 9302, 2486, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9542, 0, 3, 9302,
                                                                       2351, 9308, 2495, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9560, 0, 3, 9314,
                                                                       2363, 9320, 2522, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9578, 0, 3, 9320,
                                                                       2366, 9326, 2531, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9596, 0, 3, 9326,
                                                                       2369, 9332, 2540, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9614, 0, 3, 9332,
                                                                       2372, 9338, 2549, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9632, 0, 3, 9338,
                                                                       2375, 9344, 2558, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9650, 0, 3, 9344,
                                                                       2378, 9350, 2567, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9668, 0, 3, 9350,
                                                                       2381, 9356, 2576, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9686, 0, 3, 9356,
                                                                       2384, 9362, 2585, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9704, 0, 3, 9362,
                                                                       2387, 9368, 2594, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9722, 0, 3, 9368,
                                                                       2390, 9374, 2603, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9740, 0, 3, 9380,
                                                                       2414, 9398, 106, 112,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9776, 0, 3, 9398,
                                                                       2423, 9416, 112, 118,
                                                                       2666, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9812, 0, 3, 9416,
                                                                       2432, 9434, 118, 124,
                                                                       2684, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9848, 0, 3, 9434,
                                                                       2441, 9452, 124, 130,
                                                                       2702, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9884, 0, 3, 9452,
                                                                       2450, 9470, 130, 136,
                                                                       2720, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9920, 0, 3, 9470,
                                                                       2459, 9488, 136, 142,
                                                                       2738, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9956, 0, 3, 9488,
                                                                       2468, 9506, 142, 148,
                                                                       2756, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9992, 0, 3, 9506,
                                                                       2477, 9524, 148, 154,
                                                                       2774, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10028, 0, 3, 9524,
                                                                       2486, 9542, 154, 160,
                                                                       2792, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10064, 0, 3, 9560,
                                                                       2522, 9578, 172, 178,
                                                                       2846, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10100, 0, 3, 9578,
                                                                       2531, 9596, 178, 184,
                                                                       2864, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10136, 0, 3, 9596,
                                                                       2540, 9614, 184, 190,
                                                                       2882, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10172, 0, 3, 9614,
                                                                       2549, 9632, 190, 196,
                                                                       2900, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10208, 0, 3, 9632,
                                                                       2558, 9650, 196, 202,
                                                                       2918, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10244, 0, 3, 9650,
                                                                       2567, 9668, 202, 208,
                                                                       2936, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10280, 0, 3, 9668,
                                                                       2576, 9686, 208, 214,
                                                                       2954, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10316, 0, 3, 9686,
                                                                       2585, 9704, 214, 220,
                                                                       2972, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10352, 0, 3, 9704,
                                                                       2594, 9722, 220, 226,
                                                                       2990, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10388, 0, 3, 9740,
                                                                       2648, 9776, 238, 248,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10448, 0, 3, 9776,
                                                                       2666, 9812, 248, 258,
                                                                       3098, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10508, 0, 3, 9812,
                                                                       2684, 9848, 258, 268,
                                                                       3128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10568, 0, 3, 9848,
                                                                       2702, 9884, 268, 278,
                                                                       3158, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10628, 0, 3, 9884,
                                                                       2720, 9920, 278, 288,
                                                                       3188, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10688, 0, 3, 9920,
                                                                       2738, 9956, 288, 298,
                                                                       3218, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10748, 0, 3, 9956,
                                                                       2756, 9992, 298, 308,
                                                                       3248, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10808, 0, 3, 9992,
                                                                       2774, 10028, 308, 318,
                                                                       3278, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10868, 0, 3,
                                                                       10064, 2846, 10100, 338,
                                                                       348, 3368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10928, 0, 3,
                                                                       10100, 2864, 10136, 348,
                                                                       358, 3398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10988, 0, 3,
                                                                       10136, 2882, 10172, 358,
                                                                       368, 3428, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11048, 0, 3,
                                                                       10172, 2900, 10208, 368,
                                                                       378, 3458, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11108, 0, 3,
                                                                       10208, 2918, 10244, 378,
                                                                       388, 3488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11168, 0, 3,
                                                                       10244, 2936, 10280, 388,
                                                                       398, 3518, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11228, 0, 3,
                                                                       10280, 2954, 10316, 398,
                                                                       408, 3548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11288, 0, 3,
                                                                       10316, 2972, 10352, 408,
                                                                       418, 3578, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11348, 0, 3,
                                                                       10388, 3068, 10448, 438,
                                                                       453, 3698, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11438, 0, 3,
                                                                       10448, 3098, 10508, 453,
                                                                       468, 3743, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11528, 0, 3,
                                                                       10508, 3128, 10568, 468,
                                                                       483, 3788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11618, 0, 3,
                                                                       10568, 3158, 10628, 483,
                                                                       498, 3833, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11708, 0, 3,
                                                                       10628, 3188, 10688, 498,
                                                                       513, 3878, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11798, 0, 3,
                                                                       10688, 3218, 10748, 513,
                                                                       528, 3923, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11888, 0, 3,
                                                                       10748, 3248, 10808, 528,
                                                                       543, 3968, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11978, 0, 3,
                                                                       10868, 3368, 10928, 573,
                                                                       588, 4103, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12068, 0, 3,
                                                                       10928, 3398, 10988, 588,
                                                                       603, 4148, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12158, 0, 3,
                                                                       10988, 3428, 11048, 603,
                                                                       618, 4193, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12248, 0, 3,
                                                                       11048, 3458, 11108, 618,
                                                                       633, 4238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12338, 0, 3,
                                                                       11108, 3488, 11168, 633,
                                                                       648, 4283, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12428, 0, 3,
                                                                       11168, 3518, 11228, 648,
                                                                       663, 4328, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12518, 0, 3,
                                                                       11228, 3548, 11288, 663,
                                                                       678, 4373, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12608, 0, 3,
                                                                       11348, 3698, 11438, 708,
                                                                       729, 4544, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12734, 0, 3,
                                                                       11438, 3743, 11528, 729,
                                                                       750, 4607, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12860, 0, 3,
                                                                       11528, 3788, 11618, 750,
                                                                       771, 4670, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12986, 0, 3,
                                                                       11618, 3833, 11708, 771,
                                                                       792, 4733, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13112, 0, 3,
                                                                       11708, 3878, 11798, 792,
                                                                       813, 4796, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13238, 0, 3,
                                                                       11798, 3923, 11888, 813,
                                                                       834, 4859, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13364, 0, 3,
                                                                       11978, 4103, 12068, 876,
                                                                       897, 5048, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13490, 0, 3,
                                                                       12068, 4148, 12158, 897,
                                                                       918, 5111, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13616, 0, 3,
                                                                       12158, 4193, 12248, 918,
                                                                       939, 5174, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13742, 0, 3,
                                                                       12248, 4238, 12338, 939,
                                                                       960, 5237, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13868, 0, 3,
                                                                       12338, 4283, 12428, 960,
                                                                       981, 5300, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13994, 0, 3,
                                                                       12428, 4328, 12518, 981,
                                                                       1002, 5363, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14120, 0, 3,
                                                                       12608, 4544, 12734, 1044,
                                                                       1072, 5594, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14288, 0, 3,
                                                                       12734, 4607, 12860, 1072,
                                                                       1100, 5678, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14456, 0, 3,
                                                                       12860, 4670, 12986, 1100,
                                                                       1128, 5762, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14624, 0, 3,
                                                                       12986, 4733, 13112, 1128,
                                                                       1156, 5846, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14792, 0, 3,
                                                                       13112, 4796, 13238, 1156,
                                                                       1184, 5930, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14960, 0, 3,
                                                                       13364, 5048, 13490, 1240,
                                                                       1268, 6182, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15128, 0, 3,
                                                                       13490, 5111, 13616, 1268,
                                                                       1296, 6266, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15296, 0, 3,
                                                                       13616, 5174, 13742, 1296,
                                                                       1324, 6350, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15464, 0, 3,
                                                                       13742, 5237, 13868, 1324,
                                                                       1352, 6434, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15632, 0, 3,
                                                                       13868, 5300, 13994, 1352,
                                                                       1380, 6518, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15800, 0, 3,
                                                                       14120, 5594, 14288, 1436,
                                                                       1472, 6818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16016, 0, 3,
                                                                       14288, 5678, 14456, 1472,
                                                                       1508, 6926, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16232, 0, 3,
                                                                       14456, 5762, 14624, 1508,
                                                                       1544, 7034, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16448, 0, 3,
                                                                       14624, 5846, 14792, 1544,
                                                                       1580, 7142, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16664, 0, 3,
                                                                       14960, 6182, 15128, 1652,
                                                                       1688, 7466, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16880, 0, 3,
                                                                       15128, 6266, 15296, 1688,
                                                                       1724, 7574, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17096, 0, 3,
                                                                       15296, 6350, 15464, 1724,
                                                                       1760, 7682, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 17312, 0, 3,
                                                                       15464, 6434, 15632, 1760,
                                                                       1796, 7790, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17528, 0, 3,
                                                                       15800, 6818, 16016, 1868,
                                                                       1913, 8168, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17798, 0, 3,
                                                                       16016, 6926, 16232, 1913,
                                                                       1958, 8303, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18068, 0, 3,
                                                                       16232, 7034, 16448, 1958,
                                                                       2003, 8438, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18338, 0, 3,
                                                                       16664, 7466, 16880, 2093,
                                                                       2138, 8843, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18608, 0, 3,
                                                                       16880, 7574, 17096, 2138,
                                                                       2183, 8978, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18878, 0, 3,
                                                                       17096, 7682, 17312, 2183,
                                                                       2228, 9113, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19148, 3, 2318,
                                                                       2321, 9248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19158, 3, 2321,
                                                                       2324, 9254, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19168, 3, 2324,
                                                                       2327, 9260, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19178, 3, 2327,
                                                                       2330, 9266, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19188, 3, 2330,
                                                                       2333, 9272, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19198, 3, 2333,
                                                                       2336, 9278, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19208, 3, 2336,
                                                                       2339, 9284, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19218, 3, 2339,
                                                                       2342, 9290, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19228, 3, 2342,
                                                                       2345, 9296, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19238, 3, 2345,
                                                                       2348, 9302, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19248, 3, 2348,
                                                                       2351, 9308, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19258, 3, 2357,
                                                                       2360, 9314, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19268, 3, 2360,
                                                                       2363, 9320, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19278, 3, 2363,
                                                                       2366, 9326, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19288, 3, 2366,
                                                                       2369, 9332, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19298, 3, 2369,
                                                                       2372, 9338, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19308, 3, 2372,
                                                                       2375, 9344, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19318, 3, 2375,
                                                                       2378, 9350, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19328, 3, 2378,
                                                                       2381, 9356, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19338, 3, 2381,
                                                                       2384, 9362, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19348, 3, 2384,
                                                                       2387, 9368, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 19358, 3, 2387,
                                                                       2390, 9374, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19368, 0, 3,
                                                                       19148, 9248, 19158, 2396,
                                                                       2405, 9380, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19398, 0, 3,
                                                                       19158, 9254, 19168, 2405,
                                                                       2414, 9398, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19428, 0, 3,
                                                                       19168, 9260, 19178, 2414,
                                                                       2423, 9416, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19458, 0, 3,
                                                                       19178, 9266, 19188, 2423,
                                                                       2432, 9434, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19488, 0, 3,
                                                                       19188, 9272, 19198, 2432,
                                                                       2441, 9452, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19518, 0, 3,
                                                                       19198, 9278, 19208, 2441,
                                                                       2450, 9470, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19548, 0, 3,
                                                                       19208, 9284, 19218, 2450,
                                                                       2459, 9488, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19578, 0, 3,
                                                                       19218, 9290, 19228, 2459,
                                                                       2468, 9506, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19608, 0, 3,
                                                                       19228, 9296, 19238, 2468,
                                                                       2477, 9524, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19638, 0, 3,
                                                                       19238, 9302, 19248, 2477,
                                                                       2486, 9542, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19668, 0, 3,
                                                                       19258, 9314, 19268, 2504,
                                                                       2513, 9560, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19698, 0, 3,
                                                                       19268, 9320, 19278, 2513,
                                                                       2522, 9578, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19728, 0, 3,
                                                                       19278, 9326, 19288, 2522,
                                                                       2531, 9596, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19758, 0, 3,
                                                                       19288, 9332, 19298, 2531,
                                                                       2540, 9614, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19788, 0, 3,
                                                                       19298, 9338, 19308, 2540,
                                                                       2549, 9632, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19818, 0, 3,
                                                                       19308, 9344, 19318, 2549,
                                                                       2558, 9650, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19848, 0, 3,
                                                                       19318, 9350, 19328, 2558,
                                                                       2567, 9668, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19878, 0, 3,
                                                                       19328, 9356, 19338, 2567,
                                                                       2576, 9686, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19908, 0, 3,
                                                                       19338, 9362, 19348, 2576,
                                                                       2585, 9704, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19938, 0, 3,
                                                                       19348, 9368, 19358, 2585,
                                                                       2594, 9722, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19968, 0, 3,
                                                                       19368, 9380, 19398, 2612,
                                                                       2630, 9740, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20028, 0, 3,
                                                                       19398, 9398, 19428, 2630,
                                                                       2648, 9776, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20088, 0, 3,
                                                                       19428, 9416, 19458, 2648,
                                                                       2666, 9812, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20148, 0, 3,
                                                                       19458, 9434, 19488, 2666,
                                                                       2684, 9848, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20208, 0, 3,
                                                                       19488, 9452, 19518, 2684,
                                                                       2702, 9884, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20268, 0, 3,
                                                                       19518, 9470, 19548, 2702,
                                                                       2720, 9920, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20328, 0, 3,
                                                                       19548, 9488, 19578, 2720,
                                                                       2738, 9956, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20388, 0, 3,
                                                                       19578, 9506, 19608, 2738,
                                                                       2756, 9992, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20448, 0, 3,
                                                                       19608, 9524, 19638, 2756,
                                                                       2774, 10028, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20508, 0, 3,
                                                                       19668, 9560, 19698, 2810,
                                                                       2828, 10064, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20568, 0, 3,
                                                                       19698, 9578, 19728, 2828,
                                                                       2846, 10100, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20628, 0, 3,
                                                                       19728, 9596, 19758, 2846,
                                                                       2864, 10136, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20688, 0, 3,
                                                                       19758, 9614, 19788, 2864,
                                                                       2882, 10172, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20748, 0, 3,
                                                                       19788, 9632, 19818, 2882,
                                                                       2900, 10208, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20808, 0, 3,
                                                                       19818, 9650, 19848, 2900,
                                                                       2918, 10244, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20868, 0, 3,
                                                                       19848, 9668, 19878, 2918,
                                                                       2936, 10280, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20928, 0, 3,
                                                                       19878, 9686, 19908, 2936,
                                                                       2954, 10316, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20988, 0, 3,
                                                                       19908, 9704, 19938, 2954,
                                                                       2972, 10352, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21048, 0, 3,
                                                                       19968, 9740, 20028, 3008,
                                                                       3038, 10388, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21148, 0, 3,
                                                                       20028, 9776, 20088, 3038,
                                                                       3068, 10448, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21248, 0, 3,
                                                                       20088, 9812, 20148, 3068,
                                                                       3098, 10508, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21348, 0, 3,
                                                                       20148, 9848, 20208, 3098,
                                                                       3128, 10568, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21448, 0, 3,
                                                                       20208, 9884, 20268, 3128,
                                                                       3158, 10628, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21548, 0, 3,
                                                                       20268, 9920, 20328, 3158,
                                                                       3188, 10688, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21648, 0, 3,
                                                                       20328, 9956, 20388, 3188,
                                                                       3218, 10748, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21748, 0, 3,
                                                                       20388, 9992, 20448, 3218,
                                                                       3248, 10808, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21848, 0, 3,
                                                                       20508, 10064, 20568, 3308,
                                                                       3338, 10868, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 21948, 0, 3,
                                                                       20568, 10100, 20628, 3338,
                                                                       3368, 10928, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22048, 0, 3,
                                                                       20628, 10136, 20688, 3368,
                                                                       3398, 10988, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22148, 0, 3,
                                                                       20688, 10172, 20748, 3398,
                                                                       3428, 11048, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22248, 0, 3,
                                                                       20748, 10208, 20808, 3428,
                                                                       3458, 11108, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22348, 0, 3,
                                                                       20808, 10244, 20868, 3458,
                                                                       3488, 11168, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22448, 0, 3,
                                                                       20868, 10280, 20928, 3488,
                                                                       3518, 11228, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22548, 0, 3,
                                                                       20928, 10316, 20988, 3518,
                                                                       3548, 11288, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22648, 0, 3,
                                                                       21048, 10388, 21148, 3608,
                                                                       3653, 11348, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22798, 0, 3,
                                                                       21148, 10448, 21248, 3653,
                                                                       3698, 11438, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22948, 0, 3,
                                                                       21248, 10508, 21348, 3698,
                                                                       3743, 11528, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23098, 0, 3,
                                                                       21348, 10568, 21448, 3743,
                                                                       3788, 11618, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23248, 0, 3,
                                                                       21448, 10628, 21548, 3788,
                                                                       3833, 11708, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23398, 0, 3,
                                                                       21548, 10688, 21648, 3833,
                                                                       3878, 11798, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23548, 0, 3,
                                                                       21648, 10748, 21748, 3878,
                                                                       3923, 11888, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23698, 0, 3,
                                                                       21848, 10868, 21948, 4013,
                                                                       4058, 11978, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23848, 0, 3,
                                                                       21948, 10928, 22048, 4058,
                                                                       4103, 12068, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23998, 0, 3,
                                                                       22048, 10988, 22148, 4103,
                                                                       4148, 12158, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24148, 0, 3,
                                                                       22148, 11048, 22248, 4148,
                                                                       4193, 12248, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24298, 0, 3,
                                                                       22248, 11108, 22348, 4193,
                                                                       4238, 12338, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24448, 0, 3,
                                                                       22348, 11168, 22448, 4238,
                                                                       4283, 12428, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24598, 0, 3,
                                                                       22448, 11228, 22548, 4283,
                                                                       4328, 12518, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 24748, 0, 3,
                                                                       22648, 11348, 22798, 4418,
                                                                       4481, 12608, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 24958, 0, 3,
                                                                       22798, 11438, 22948, 4481,
                                                                       4544, 12734, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25168, 0, 3,
                                                                       22948, 11528, 23098, 4544,
                                                                       4607, 12860, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25378, 0, 3,
                                                                       23098, 11618, 23248, 4607,
                                                                       4670, 12986, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25588, 0, 3,
                                                                       23248, 11708, 23398, 4670,
                                                                       4733, 13112, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25798, 0, 3,
                                                                       23398, 11798, 23548, 4733,
                                                                       4796, 13238, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26008, 0, 3,
                                                                       23698, 11978, 23848, 4922,
                                                                       4985, 13364, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26218, 0, 3,
                                                                       23848, 12068, 23998, 4985,
                                                                       5048, 13490, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26428, 0, 3,
                                                                       23998, 12158, 24148, 5048,
                                                                       5111, 13616, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26638, 0, 3,
                                                                       24148, 12248, 24298, 5111,
                                                                       5174, 13742, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26848, 0, 3,
                                                                       24298, 12338, 24448, 5174,
                                                                       5237, 13868, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 27058, 0, 3,
                                                                       24448, 12428, 24598, 5237,
                                                                       5300, 13994, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27268, 0, 3,
                                                                       24748, 12608, 24958, 5426,
                                                                       5510, 14120, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27548, 0, 3,
                                                                       24958, 12734, 25168, 5510,
                                                                       5594, 14288, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27828, 0, 3,
                                                                       25168, 12860, 25378, 5594,
                                                                       5678, 14456, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28108, 0, 3,
                                                                       25378, 12986, 25588, 5678,
                                                                       5762, 14624, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28388, 0, 3,
                                                                       25588, 13112, 25798, 5762,
                                                                       5846, 14792, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28668, 0, 3,
                                                                       26008, 13364, 26218, 6014,
                                                                       6098, 14960, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28948, 0, 3,
                                                                       26218, 13490, 26428, 6098,
                                                                       6182, 15128, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 29228, 0, 3,
                                                                       26428, 13616, 26638, 6182,
                                                                       6266, 15296, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 29508, 0, 3,
                                                                       26638, 13742, 26848, 6266,
                                                                       6350, 15464, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 29788, 0, 3,
                                                                       26848, 13868, 27058, 6350,
                                                                       6434, 15632, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 30068, 0, 3,
                                                                       27268, 14120, 27548, 6602,
                                                                       6710, 15800, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 30428, 0, 3,
                                                                       27548, 14288, 27828, 6710,
                                                                       6818, 16016, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 30788, 0, 3,
                                                                       27828, 14456, 28108, 6818,
                                                                       6926, 16232, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 31148, 0, 3,
                                                                       28108, 14624, 28388, 6926,
                                                                       7034, 16448, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 31508, 0, 3,
                                                                       28668, 14960, 28948, 7250,
                                                                       7358, 16664, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 31868, 0, 3,
                                                                       28948, 15128, 29228, 7358,
                                                                       7466, 16880, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32228, 0, 3,
                                                                       29228, 15296, 29508, 7466,
                                                                       7574, 17096, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 32588, 0, 3,
                                                                       29508, 15464, 29788, 7574,
                                                                       7682, 17312, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 32948, 0, 3,
                                                                       30068, 15800, 30428, 7898,
                                                                       8033, 17528, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 33398, 0, 3,
                                                                       30428, 16016, 30788, 8033,
                                                                       8168, 17798, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 33848, 0, 3,
                                                                       30788, 16232, 31148, 8168,
                                                                       8303, 18068, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 34298, 0, 3,
                                                                       31508, 16664, 31868, 8573,
                                                                       8708, 18338, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 34748, 0, 3,
                                                                       31868, 16880, 32228, 8708,
                                                                       8843, 18608, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 35198, 0, 3,
                                                                       32228, 17096, 32588, 8843,
                                                                       8978, 18878, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35648, 3, 9248,
                                                                       9254, 19168, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35663, 3, 9254,
                                                                       9260, 19178, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35678, 3, 9260,
                                                                       9266, 19188, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35693, 3, 9266,
                                                                       9272, 19198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35708, 3, 9272,
                                                                       9278, 19208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35723, 3, 9278,
                                                                       9284, 19218, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35738, 3, 9284,
                                                                       9290, 19228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35753, 3, 9290,
                                                                       9296, 19238, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35768, 3, 9296,
                                                                       9302, 19248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35783, 3, 9314,
                                                                       9320, 19278, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35798, 3, 9320,
                                                                       9326, 19288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35813, 3, 9326,
                                                                       9332, 19298, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35828, 3, 9332,
                                                                       9338, 19308, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35843, 3, 9338,
                                                                       9344, 19318, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35858, 3, 9344,
                                                                       9350, 19328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35873, 3, 9350,
                                                                       9356, 19338, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35888, 3, 9356,
                                                                       9362, 19348, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35903, 3, 9362,
                                                                       9368, 19358, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35918, 0, 3,
                                                                       35648, 19168, 35663, 9380,
                                                                       9398, 19428, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35963, 0, 3,
                                                                       35663, 19178, 35678, 9398,
                                                                       9416, 19458, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36008, 0, 3,
                                                                       35678, 19188, 35693, 9416,
                                                                       9434, 19488, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36053, 0, 3,
                                                                       35693, 19198, 35708, 9434,
                                                                       9452, 19518, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36098, 0, 3,
                                                                       35708, 19208, 35723, 9452,
                                                                       9470, 19548, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36143, 0, 3,
                                                                       35723, 19218, 35738, 9470,
                                                                       9488, 19578, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36188, 0, 3,
                                                                       35738, 19228, 35753, 9488,
                                                                       9506, 19608, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36233, 0, 3,
                                                                       35753, 19238, 35768, 9506,
                                                                       9524, 19638, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36278, 0, 3,
                                                                       35783, 19278, 35798, 9560,
                                                                       9578, 19728, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36323, 0, 3,
                                                                       35798, 19288, 35813, 9578,
                                                                       9596, 19758, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36368, 0, 3,
                                                                       35813, 19298, 35828, 9596,
                                                                       9614, 19788, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36413, 0, 3,
                                                                       35828, 19308, 35843, 9614,
                                                                       9632, 19818, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36458, 0, 3,
                                                                       35843, 19318, 35858, 9632,
                                                                       9650, 19848, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36503, 0, 3,
                                                                       35858, 19328, 35873, 9650,
                                                                       9668, 19878, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36548, 0, 3,
                                                                       35873, 19338, 35888, 9668,
                                                                       9686, 19908, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36593, 0, 3,
                                                                       35888, 19348, 35903, 9686,
                                                                       9704, 19938, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36638, 0, 3,
                                                                       35918, 19428, 35963, 9740,
                                                                       9776, 20088, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36728, 0, 3,
                                                                       35963, 19458, 36008, 9776,
                                                                       9812, 20148, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36818, 0, 3,
                                                                       36008, 19488, 36053, 9812,
                                                                       9848, 20208, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36908, 0, 3,
                                                                       36053, 19518, 36098, 9848,
                                                                       9884, 20268, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 36998, 0, 3,
                                                                       36098, 19548, 36143, 9884,
                                                                       9920, 20328, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37088, 0, 3,
                                                                       36143, 19578, 36188, 9920,
                                                                       9956, 20388, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37178, 0, 3,
                                                                       36188, 19608, 36233, 9956,
                                                                       9992, 20448, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37268, 0, 3,
                                                                       36278, 19728, 36323,
                                                                       10064, 10100, 20628,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37358, 0, 3,
                                                                       36323, 19758, 36368,
                                                                       10100, 10136, 20688,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37448, 0, 3,
                                                                       36368, 19788, 36413,
                                                                       10136, 10172, 20748,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37538, 0, 3,
                                                                       36413, 19818, 36458,
                                                                       10172, 10208, 20808,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37628, 0, 3,
                                                                       36458, 19848, 36503,
                                                                       10208, 10244, 20868,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37718, 0, 3,
                                                                       36503, 19878, 36548,
                                                                       10244, 10280, 20928,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37808, 0, 3,
                                                                       36548, 19908, 36593,
                                                                       10280, 10316, 20988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37898, 0, 3,
                                                                       36638, 20088, 36728,
                                                                       10388, 10448, 21248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38048, 0, 3,
                                                                       36728, 20148, 36818,
                                                                       10448, 10508, 21348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38198, 0, 3,
                                                                       36818, 20208, 36908,
                                                                       10508, 10568, 21448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38348, 0, 3,
                                                                       36908, 20268, 36998,
                                                                       10568, 10628, 21548,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38498, 0, 3,
                                                                       36998, 20328, 37088,
                                                                       10628, 10688, 21648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38648, 0, 3,
                                                                       37088, 20388, 37178,
                                                                       10688, 10748, 21748,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38798, 0, 3,
                                                                       37268, 20628, 37358,
                                                                       10868, 10928, 22048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38948, 0, 3,
                                                                       37358, 20688, 37448,
                                                                       10928, 10988, 22148,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 39098, 0, 3,
                                                                       37448, 20748, 37538,
                                                                       10988, 11048, 22248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 39248, 0, 3,
                                                                       37538, 20808, 37628,
                                                                       11048, 11108, 22348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 39398, 0, 3,
                                                                       37628, 20868, 37718,
                                                                       11108, 11168, 22448,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 39548, 0, 3,
                                                                       37718, 20928, 37808,
                                                                       11168, 11228, 22548,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 39698, 0, 3,
                                                                       37898, 21248, 38048,
                                                                       11348, 11438, 22948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 39923, 0, 3,
                                                                       38048, 21348, 38198,
                                                                       11438, 11528, 23098,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 40148, 0, 3,
                                                                       38198, 21448, 38348,
                                                                       11528, 11618, 23248,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 40373, 0, 3,
                                                                       38348, 21548, 38498,
                                                                       11618, 11708, 23398,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 40598, 0, 3,
                                                                       38498, 21648, 38648,
                                                                       11708, 11798, 23548,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 40823, 0, 3,
                                                                       38798, 22048, 38948,
                                                                       11978, 12068, 23998,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 41048, 0, 3,
                                                                       38948, 22148, 39098,
                                                                       12068, 12158, 24148,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 41273, 0, 3,
                                                                       39098, 22248, 39248,
                                                                       12158, 12248, 24298,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 41498, 0, 3,
                                                                       39248, 22348, 39398,
                                                                       12248, 12338, 24448,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 41723, 0, 3,
                                                                       39398, 22448, 39548,
                                                                       12338, 12428, 24598,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 41948, 0, 3,
                                                                       39698, 22948, 39923,
                                                                       12608, 12734, 25168,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 42263, 0, 3,
                                                                       39923, 23098, 40148,
                                                                       12734, 12860, 25378,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 42578, 0, 3,
                                                                       40148, 23248, 40373,
                                                                       12860, 12986, 25588,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 42893, 0, 3,
                                                                       40373, 23398, 40598,
                                                                       12986, 13112, 25798,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 43208, 0, 3,
                                                                       40823, 23998, 41048,
                                                                       13364, 13490, 26428,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 43523, 0, 3,
                                                                       41048, 24148, 41273,
                                                                       13490, 13616, 26638,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 43838, 0, 3,
                                                                       41273, 24298, 41498,
                                                                       13616, 13742, 26848,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 44153, 0, 3,
                                                                       41498, 24448, 41723,
                                                                       13742, 13868, 27058,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 44468, 0, 3,
                                                                       41948, 25168, 42263,
                                                                       14120, 14288, 27828,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 44888, 0, 3,
                                                                       42263, 25378, 42578,
                                                                       14288, 14456, 28108,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 45308, 0, 3,
                                                                       42578, 25588, 42893,
                                                                       14456, 14624, 28388,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 45728, 0, 3,
                                                                       43208, 26428, 43523,
                                                                       14960, 15128, 29228,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 46148, 0, 3,
                                                                       43523, 26638, 43838,
                                                                       15128, 15296, 29508,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 46568, 0, 3,
                                                                       43838, 26848, 44153,
                                                                       15296, 15464, 29788,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 46988, 0, 3,
                                                                       44468, 27828, 44888,
                                                                       15800, 16016, 30788,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 47528, 0, 3,
                                                                       44888, 28108, 45308,
                                                                       16016, 16232, 31148,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 48068, 0, 3,
                                                                       45728, 29228, 46148,
                                                                       16664, 16880, 32228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 48608, 0, 3,
                                                                       46148, 29508, 46568,
                                                                       16880, 17096, 32588,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 49148, 0, 3,
                                                                       46988, 30788, 47528,
                                                                       17528, 17798, 33848,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 49823, 0, 3,
                                                                       48068, 32228, 48608,
                                                                       18338, 18608, 35198,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50498, 3, 19148,
                                                                       19158, 35648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50519, 3, 19158,
                                                                       19168, 35663, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50540, 3, 19168,
                                                                       19178, 35678, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50561, 3, 19178,
                                                                       19188, 35693, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50582, 3, 19188,
                                                                       19198, 35708, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50603, 3, 19198,
                                                                       19208, 35723, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50624, 3, 19208,
                                                                       19218, 35738, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50645, 3, 19218,
                                                                       19228, 35753, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50666, 3, 19228,
                                                                       19238, 35768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50687, 3, 19258,
                                                                       19268, 35783, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50708, 3, 19268,
                                                                       19278, 35798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50729, 3, 19278,
                                                                       19288, 35813, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50750, 3, 19288,
                                                                       19298, 35828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50771, 3, 19298,
                                                                       19308, 35843, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50792, 3, 19308,
                                                                       19318, 35858, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50813, 3, 19318,
                                                                       19328, 35873, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50834, 3, 19328,
                                                                       19338, 35888, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50855, 3, 19338,
                                                                       19348, 35903, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 50876, 0, 3,
                                                                       50498, 35648, 50519,
                                                                       19368, 19398, 35918,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 50939, 0, 3,
                                                                       50519, 35663, 50540,
                                                                       19398, 19428, 35963,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51002, 0, 3,
                                                                       50540, 35678, 50561,
                                                                       19428, 19458, 36008,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51065, 0, 3,
                                                                       50561, 35693, 50582,
                                                                       19458, 19488, 36053,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51128, 0, 3,
                                                                       50582, 35708, 50603,
                                                                       19488, 19518, 36098,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51191, 0, 3,
                                                                       50603, 35723, 50624,
                                                                       19518, 19548, 36143,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51254, 0, 3,
                                                                       50624, 35738, 50645,
                                                                       19548, 19578, 36188,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51317, 0, 3,
                                                                       50645, 35753, 50666,
                                                                       19578, 19608, 36233,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51380, 0, 3,
                                                                       50687, 35783, 50708,
                                                                       19668, 19698, 36278,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51443, 0, 3,
                                                                       50708, 35798, 50729,
                                                                       19698, 19728, 36323,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51506, 0, 3,
                                                                       50729, 35813, 50750,
                                                                       19728, 19758, 36368,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51569, 0, 3,
                                                                       50750, 35828, 50771,
                                                                       19758, 19788, 36413,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51632, 0, 3,
                                                                       50771, 35843, 50792,
                                                                       19788, 19818, 36458,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51695, 0, 3,
                                                                       50792, 35858, 50813,
                                                                       19818, 19848, 36503,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51758, 0, 3,
                                                                       50813, 35873, 50834,
                                                                       19848, 19878, 36548,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 51821, 0, 3,
                                                                       50834, 35888, 50855,
                                                                       19878, 19908, 36593,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 51884, 0, 3,
                                                                       50876, 35918, 50939,
                                                                       19968, 20028, 36638,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52010, 0, 3,
                                                                       50939, 35963, 51002,
                                                                       20028, 20088, 36728,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52136, 0, 3,
                                                                       51002, 36008, 51065,
                                                                       20088, 20148, 36818,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52262, 0, 3,
                                                                       51065, 36053, 51128,
                                                                       20148, 20208, 36908,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52388, 0, 3,
                                                                       51128, 36098, 51191,
                                                                       20208, 20268, 36998,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52514, 0, 3,
                                                                       51191, 36143, 51254,
                                                                       20268, 20328, 37088,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52640, 0, 3,
                                                                       51254, 36188, 51317,
                                                                       20328, 20388, 37178,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52766, 0, 3,
                                                                       51380, 36278, 51443,
                                                                       20508, 20568, 37268,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 52892, 0, 3,
                                                                       51443, 36323, 51506,
                                                                       20568, 20628, 37358,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53018, 0, 3,
                                                                       51506, 36368, 51569,
                                                                       20628, 20688, 37448,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53144, 0, 3,
                                                                       51569, 36413, 51632,
                                                                       20688, 20748, 37538,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53270, 0, 3,
                                                                       51632, 36458, 51695,
                                                                       20748, 20808, 37628,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53396, 0, 3,
                                                                       51695, 36503, 51758,
                                                                       20808, 20868, 37718,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 53522, 0, 3,
                                                                       51758, 36548, 51821,
                                                                       20868, 20928, 37808,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 53648, 0, 3,
                                                                       51884, 36638, 52010,
                                                                       21048, 21148, 37898,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 53858, 0, 3,
                                                                       52010, 36728, 52136,
                                                                       21148, 21248, 38048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54068, 0, 3,
                                                                       52136, 36818, 52262,
                                                                       21248, 21348, 38198,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54278, 0, 3,
                                                                       52262, 36908, 52388,
                                                                       21348, 21448, 38348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54488, 0, 3,
                                                                       52388, 36998, 52514,
                                                                       21448, 21548, 38498,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54698, 0, 3,
                                                                       52514, 37088, 52640,
                                                                       21548, 21648, 38648,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54908, 0, 3,
                                                                       52766, 37268, 52892,
                                                                       21848, 21948, 38798,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 55118, 0, 3,
                                                                       52892, 37358, 53018,
                                                                       21948, 22048, 38948,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 55328, 0, 3,
                                                                       53018, 37448, 53144,
                                                                       22048, 22148, 39098,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 55538, 0, 3,
                                                                       53144, 37538, 53270,
                                                                       22148, 22248, 39248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 55748, 0, 3,
                                                                       53270, 37628, 53396,
                                                                       22248, 22348, 39398,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 55958, 0, 3,
                                                                       53396, 37718, 53522,
                                                                       22348, 22448, 39548,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 56168, 0, 3,
                                                                       53648, 37898, 53858,
                                                                       22648, 22798, 39698,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 56483, 0, 3,
                                                                       53858, 38048, 54068,
                                                                       22798, 22948, 39923,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 56798, 0, 3,
                                                                       54068, 38198, 54278,
                                                                       22948, 23098, 40148,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 57113, 0, 3,
                                                                       54278, 38348, 54488,
                                                                       23098, 23248, 40373,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 57428, 0, 3,
                                                                       54488, 38498, 54698,
                                                                       23248, 23398, 40598,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 57743, 0, 3,
                                                                       54908, 38798, 55118,
                                                                       23698, 23848, 40823,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 58058, 0, 3,
                                                                       55118, 38948, 55328,
                                                                       23848, 23998, 41048,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 58373, 0, 3,
                                                                       55328, 39098, 55538,
                                                                       23998, 24148, 41273,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 58688, 0, 3,
                                                                       55538, 39248, 55748,
                                                                       24148, 24298, 41498,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 59003, 0, 3,
                                                                       55748, 39398, 55958,
                                                                       24298, 24448, 41723,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 59318, 0, 3,
                                                                       56168, 39698, 56483,
                                                                       24748, 24958, 41948,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 59759, 0, 3,
                                                                       56483, 39923, 56798,
                                                                       24958, 25168, 42263,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 60200, 0, 3,
                                                                       56798, 40148, 57113,
                                                                       25168, 25378, 42578,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 60641, 0, 3,
                                                                       57113, 40373, 57428,
                                                                       25378, 25588, 42893,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 61082, 0, 3,
                                                                       57743, 40823, 58058,
                                                                       26008, 26218, 43208,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 61523, 0, 3,
                                                                       58058, 41048, 58373,
                                                                       26218, 26428, 43523,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 61964, 0, 3,
                                                                       58373, 41273, 58688,
                                                                       26428, 26638, 43838,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 62405, 0, 3,
                                                                       58688, 41498, 59003,
                                                                       26638, 26848, 44153,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 62846, 0, 3,
                                                                       59318, 41948, 59759,
                                                                       27268, 27548, 44468,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 63434, 0, 3,
                                                                       59759, 42263, 60200,
                                                                       27548, 27828, 44888,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 64022, 0, 3,
                                                                       60200, 42578, 60641,
                                                                       27828, 28108, 45308,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 64610, 0, 3,
                                                                       61082, 43208, 61523,
                                                                       28668, 28948, 45728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 65198, 0, 3,
                                                                       61523, 43523, 61964,
                                                                       28948, 29228, 46148,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 65786, 0, 3,
                                                                       61964, 43838, 62405,
                                                                       29228, 29508, 46568,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 66374, 0, 3,
                                                                       62846, 44468, 63434,
                                                                       30068, 30428, 46988,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 67130, 0, 3,
                                                                       63434, 44888, 64022,
                                                                       30428, 30788, 47528,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 67886, 0, 3,
                                                                       64610, 45728, 65198,
                                                                       31508, 31868, 48068,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 68642, 0, 3,
                                                                       65198, 46148, 65786,
                                                                       31868, 32228, 48608,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 69398, 0, 3,
                                                                       66374, 46988, 67130,
                                                                       32948, 33398, 49148,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 70343, 0, 3,
                                                                       67886, 48068, 68642,
                                                                       34298, 34748, 49823,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_f_x(buffer, 71288, 51884, 56168, 1, 21, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 71498, 51884, 56168, 1, 21, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 71708, 51884, 56168, 1, 21, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 71918, 52766, 57743, 1, 21, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 72128, 52766, 57743, 1, 21, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 72338, 52766, 57743, 1, 21, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 72548, 53648, 59318, 1, 21, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 72863, 53648, 59318, 1, 21, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 73178, 53648, 59318, 1, 21, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 73493, 54908, 61082, 1, 21, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 73808, 54908, 61082, 1, 21, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 74123, 54908, 61082, 1, 21, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 74438, 56168, 62846, 1, 21, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 74879, 56168, 62846, 1, 21, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 75320, 56168, 62846, 1, 21, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 75761, 57743, 64610, 1, 21, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 76202, 57743, 64610, 1, 21, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 76643, 57743, 64610, 1, 21, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 77084, 59318, 66374, 1, 21, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 77672, 59318, 66374, 1, 21, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 78260, 59318, 66374, 1, 21, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 78848, 61082, 67886, 1, 21, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 79436, 61082, 67886, 1, 21, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 80024, 61082, 67886, 1, 21, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 80612, 62846, 69398, 1, 21, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 81368, 62846, 69398, 1, 21, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 82124, 62846, 69398, 1, 21, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 82880, 64610, 70343, 1, 21, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 83636, 64610, 70343, 1, 21, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 84392, 64610, 70343, 1, 21, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 85148, 71288, 210, ncols);

                    simdfunc::contract_primitives(buffer, 85468, 71498, 210, ncols);

                    simdfunc::contract_primitives(buffer, 85788, 71708, 210, ncols);

                    simdfunc::contract_primitives(buffer, 86108, 53648, 210, ncols);

                    simdfunc::contract_primitives(buffer, 86428, 71918, 210, ncols);

                    simdfunc::contract_primitives(buffer, 86748, 72128, 210, ncols);

                    simdfunc::contract_primitives(buffer, 87068, 72338, 210, ncols);

                    simdfunc::contract_primitives(buffer, 87388, 54908, 210, ncols);

                    simdfunc::contract_primitives(buffer, 87708, 72548, 315, ncols);

                    simdfunc::contract_primitives(buffer, 88188, 72863, 315, ncols);

                    simdfunc::contract_primitives(buffer, 88668, 73178, 315, ncols);

                    simdfunc::contract_primitives(buffer, 89148, 56168, 315, ncols);

                    simdfunc::contract_primitives(buffer, 89628, 73493, 315, ncols);

                    simdfunc::contract_primitives(buffer, 90108, 73808, 315, ncols);

                    simdfunc::contract_primitives(buffer, 90588, 74123, 315, ncols);

                    simdfunc::contract_primitives(buffer, 91068, 57743, 315, ncols);

                    simdfunc::contract_primitives(buffer, 91548, 74438, 441, ncols);

                    simdfunc::contract_primitives(buffer, 92220, 74879, 441, ncols);

                    simdfunc::contract_primitives(buffer, 92892, 75320, 441, ncols);

                    simdfunc::contract_primitives(buffer, 93564, 59318, 441, ncols);

                    simdfunc::contract_primitives(buffer, 94236, 75761, 441, ncols);

                    simdfunc::contract_primitives(buffer, 94908, 76202, 441, ncols);

                    simdfunc::contract_primitives(buffer, 95580, 76643, 441, ncols);

                    simdfunc::contract_primitives(buffer, 96252, 61082, 441, ncols);

                    simdfunc::contract_primitives(buffer, 96924, 77084, 588, ncols);

                    simdfunc::contract_primitives(buffer, 97820, 77672, 588, ncols);

                    simdfunc::contract_primitives(buffer, 98716, 78260, 588, ncols);

                    simdfunc::contract_primitives(buffer, 99612, 62846, 588, ncols);

                    simdfunc::contract_primitives(buffer, 100508, 78848, 588, ncols);

                    simdfunc::contract_primitives(buffer, 101404, 79436, 588, ncols);

                    simdfunc::contract_primitives(buffer, 102300, 80024, 588, ncols);

                    simdfunc::contract_primitives(buffer, 103196, 64610, 588, ncols);

                    simdfunc::contract_primitives(buffer, 104092, 80612, 756, ncols);

                    simdfunc::contract_primitives(buffer, 105244, 81368, 756, ncols);

                    simdfunc::contract_primitives(buffer, 106396, 82124, 756, ncols);

                    simdfunc::contract_primitives(buffer, 107548, 82880, 756, ncols);

                    simdfunc::contract_primitives(buffer, 108700, 83636, 756, ncols);

                    simdfunc::contract_primitives(buffer, 109852, 84392, 756, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 85358, 85148, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 85678, 85468, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 85998, 85788, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 86318, 86108, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 86638, 86428, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 86958, 86748, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 87278, 87068, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 87598, 87388, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 88023, 87708, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 88503, 88188, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 88983, 88668, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 89463, 89148, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 89943, 89628, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 90423, 90108, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 90903, 90588, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 91383, 91068, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 91989, 91548, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 92661, 92220, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 93333, 92892, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 94005, 93564, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 94677, 94236, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 95349, 94908, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 96021, 95580, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 96693, 96252, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 97512, 96924, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 98408, 97820, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 99304, 98716, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 100200, 99612, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 101096, 100508, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 101992, 101404, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 102888, 102300, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 103784, 103196, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 104848, 104092, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 106000, 105244, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 107152, 106396, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 108304, 107548, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 109456, 108700, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 110608, 109852, 36, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 111004, 85358, 86318,
                                                       88023, 11, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 111334, 85678, 86318,
                                                       88503, 11, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 111664, 85998, 86318,
                                                       88983, 11, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 111994, 86318, 89463, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 112324, 86638, 87598,
                                                       89943, 11, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 112654, 86958, 87598,
                                                       90423, 11, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 112984, 87278, 87598,
                                                       90903, 11, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 113314, 87598, 91383, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 113644, 88023, 89463,
                                                       91989, 11, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 114139, 88503, 89463,
                                                       92661, 11, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 114634, 88983, 89463,
                                                       93333, 11, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 115129, 89463, 94005, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 115624, 89943, 91383,
                                                       94677, 11, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 116119, 90423, 91383,
                                                       95349, 11, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 116614, 90903, 91383,
                                                       96021, 11, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 117109, 91383, 96693, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 117604, 91989, 94005,
                                                       97512, 11, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 118297, 92661, 94005,
                                                       98408, 11, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 118990, 93333, 94005,
                                                       99304, 11, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 119683, 94005, 100200, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 120376, 94677, 96693,
                                                       101096, 11, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 121069, 95349, 96693,
                                                       101992, 11, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 121762, 96021, 96693,
                                                       102888, 11, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 122455, 96693, 103784, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 123148, 97512,
                                                       100200, 104848, 11, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 124072, 98408,
                                                       100200, 106000, 11, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 124996, 99304,
                                                       100200, 107152, 11, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 125920, 101096,
                                                       103784, 108304, 11, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 126844, 101992,
                                                       103784, 109456, 11, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 127768, 102888,
                                                       103784, 110608, 11, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 128692, 111004,
                                                       111994, 113644, 11, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 129352, 111334,
                                                       111994, 114139, 11, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 130012, 111664,
                                                       111994, 114634, 11, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 130672, 111994, 115129, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 131332, 112324,
                                                       113314, 115624, 11, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 131992, 112654,
                                                       113314, 116119, 11, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 132652, 112984,
                                                       113314, 116614, 11, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 133312, 113314, 117109, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 133972, 113644,
                                                       115129, 117604, 11, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 134962, 114139,
                                                       115129, 118297, 11, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 135952, 114634,
                                                       115129, 118990, 11, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 136942, 115129, 119683, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 137932, 115624,
                                                       117109, 120376, 11, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 138922, 116119,
                                                       117109, 121069, 11, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 139912, 116614,
                                                       117109, 121762, 11, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 140902, 117109, 122455, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 141892, 117604,
                                                       119683, 123148, 11, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 143278, 118297,
                                                       119683, 124072, 11, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 144664, 118990,
                                                       119683, 124996, 11, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 146050, 120376,
                                                       122455, 125920, 11, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 147436, 121069,
                                                       122455, 126844, 11, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 148822, 121762,
                                                       122455, 127768, 11, nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 150208, 128692,
                                                       130672, 133972, 11, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 151308, 129352,
                                                       130672, 134962, 11, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 152408, 130012,
                                                       130672, 135952, 11, nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 153508, 130672, 136942, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 154608, 131332,
                                                       133312, 137932, 11, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 155708, 131992,
                                                       133312, 138922, 11, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 156808, 132652,
                                                       133312, 139912, 11, nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 157908, 133312, 140902, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 159008, 133972,
                                                       136942, 141892, 11, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 160658, 134962,
                                                       136942, 143278, 11, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 162308, 135952,
                                                       136942, 144664, 11, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 163958, 137932,
                                                       140902, 146050, 11, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 165608, 138922,
                                                       140902, 147436, 11, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 167258, 139912,
                                                       140902, 148822, 11, nmax);

        simdtrf::compute_hrr_geom_100x_fg_out_of_first(buffer, coordinates, 168908, 150208,
                                                       153508, 159008, 11, nmax);

        simdtrf::compute_hrr_geom_100y_fg_out_of_first(buffer, coordinates, 170558, 151308,
                                                       153508, 160658, 11, nmax);

        simdtrf::compute_hrr_geom_100z_fg_out_of_first(buffer, coordinates, 172208, 152408,
                                                       153508, 162308, 11, nmax);

        simdtrf::compute_hrr_geom_100x_fg_out_of_first(buffer, coordinates, 173858, 154608,
                                                       157908, 163958, 11, nmax);

        simdtrf::compute_hrr_geom_100y_fg_out_of_first(buffer, coordinates, 175508, 155708,
                                                       157908, 165608, 11, nmax);

        simdtrf::compute_hrr_geom_100z_fg_out_of_first(buffer, coordinates, 177158, 156808,
                                                       157908, 167258, 11, nmax);

        simdtrf::transform_g_inner(buffer, 178808, 173858, 10, 11, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 178808, 99, nmax);

        simdtrf::transform_g_inner(buffer, 178808, 175508, 10, 11, nmax);

        simdtrf::transform_f_outer(values + 693 * nvalues + n * npairs, nvalues, buffer, 178808,
                                   99, nmax);

        simdtrf::transform_g_inner(buffer, 178808, 177158, 10, 11, nmax);

        simdtrf::transform_f_outer(values + 1386 * nvalues + n * npairs, nvalues, buffer, 178808,
                                   99, nmax);

        simdtrf::transform_g_inner(buffer, 178808, 168908, 10, 11, nmax);

        simdtrf::transform_f_outer(values + 2079 * nvalues + n * npairs, nvalues, buffer, 178808,
                                   99, nmax);

        simdtrf::transform_g_inner(buffer, 178808, 170558, 10, 11, nmax);

        simdtrf::transform_f_outer(values + 2772 * nvalues + n * npairs, nvalues, buffer, 178808,
                                   99, nmax);

        simdtrf::transform_g_inner(buffer, 178808, 172208, 10, 11, nmax);

        simdtrf::transform_f_outer(values + 3465 * nvalues + n * npairs, nvalues, buffer, 178808,
                                   99, nmax);
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
