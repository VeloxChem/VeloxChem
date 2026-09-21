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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecFGI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
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
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_fgi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_fgi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 246438, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 4914 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 246438, 132008, 32660, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 14,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 22, 3, 14,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 77, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 89, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 95, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 101, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 104, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 107, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 110, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 116, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 119, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 7, 8,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 8, 9,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 134, 0, 3, 9, 10,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 140, 0, 3, 10, 11,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 146, 0, 3, 11, 12,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 12, 13,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 13, 14,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 164, 0, 3, 14, 15,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 170, 0, 3, 15, 16,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 176, 0, 3, 16, 17,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 17, 18,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 18, 19,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 194, 0, 3, 19, 20,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 200, 0, 3, 23, 24,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 206, 0, 3, 24, 25,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 25, 26,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 26, 27,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 224, 0, 3, 27, 28,
                                                                       92, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 230, 0, 3, 28, 29,
                                                                       95, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 236, 0, 3, 29, 30,
                                                                       98, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 242, 0, 3, 30, 31,
                                                                       101, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 31, 32,
                                                                       104, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 254, 0, 3, 32, 33,
                                                                       107, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 260, 0, 3, 33, 34,
                                                                       110, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 266, 0, 3, 34, 35,
                                                                       113, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 272, 0, 3, 35, 36,
                                                                       116, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 38, 41,
                                                                       122, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 41, 44,
                                                                       128, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 44, 47,
                                                                       134, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 47, 50,
                                                                       140, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 50, 53,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 53, 56,
                                                                       152, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 56, 59,
                                                                       158, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 59, 62,
                                                                       164, 170, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 62, 65,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 65, 68,
                                                                       176, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 68, 71,
                                                                       182, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 71, 74,
                                                                       188, 194, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 80, 83,
                                                                       200, 206, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 83, 86,
                                                                       206, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 86, 89,
                                                                       212, 218, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 89, 92,
                                                                       218, 224, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 92, 95,
                                                                       224, 230, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 95, 98,
                                                                       230, 236, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 98,
                                                                       101, 236, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 101,
                                                                       104, 242, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 478, 0, 3, 104,
                                                                       107, 248, 254, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 488, 0, 3, 107,
                                                                       110, 254, 260, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 110,
                                                                       113, 260, 266, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 508, 0, 3, 113,
                                                                       116, 266, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 122,
                                                                       128, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 533, 0, 3, 128,
                                                                       134, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 548, 0, 3, 134,
                                                                       140, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 563, 0, 3, 140,
                                                                       146, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 578, 0, 3, 146,
                                                                       152, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 593, 0, 3, 152,
                                                                       158, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 608, 0, 3, 158,
                                                                       164, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 623, 0, 3, 164,
                                                                       170, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 170,
                                                                       176, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 653, 0, 3, 176,
                                                                       182, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 668, 0, 3, 182,
                                                                       188, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 683, 0, 3, 200,
                                                                       206, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 698, 0, 3, 206,
                                                                       212, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 713, 0, 3, 212,
                                                                       218, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 728, 0, 3, 218,
                                                                       224, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 743, 0, 3, 224,
                                                                       230, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 758, 0, 3, 230,
                                                                       236, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 773, 0, 3, 236,
                                                                       242, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 788, 0, 3, 242,
                                                                       248, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 803, 0, 3, 248,
                                                                       254, 478, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 818, 0, 3, 254,
                                                                       260, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 833, 0, 3, 260,
                                                                       266, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 848, 0, 3, 278,
                                                                       288, 518, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 869, 0, 3, 288,
                                                                       298, 533, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 890, 0, 3, 298,
                                                                       308, 548, 563, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 911, 0, 3, 308,
                                                                       318, 563, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 932, 0, 3, 318,
                                                                       328, 578, 593, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 953, 0, 3, 328,
                                                                       338, 593, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 974, 0, 3, 338,
                                                                       348, 608, 623, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 995, 0, 3, 348,
                                                                       358, 623, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1016, 0, 3, 358,
                                                                       368, 638, 653, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1037, 0, 3, 368,
                                                                       378, 653, 668, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 398,
                                                                       408, 683, 698, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1079, 0, 3, 408,
                                                                       418, 698, 713, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 418,
                                                                       428, 713, 728, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1121, 0, 3, 428,
                                                                       438, 728, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 438,
                                                                       448, 743, 758, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1163, 0, 3, 448,
                                                                       458, 758, 773, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 458,
                                                                       468, 773, 788, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1205, 0, 3, 468,
                                                                       478, 788, 803, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1226, 0, 3, 478,
                                                                       488, 803, 818, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1247, 0, 3, 488,
                                                                       498, 818, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 518,
                                                                       533, 848, 869, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 533,
                                                                       548, 869, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 548,
                                                                       563, 890, 911, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 563,
                                                                       578, 911, 932, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 578,
                                                                       593, 932, 953, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 593,
                                                                       608, 953, 974, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 608,
                                                                       623, 974, 995, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1464, 0, 3, 623,
                                                                       638, 995, 1016, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 638,
                                                                       653, 1016, 1037, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 683,
                                                                       698, 1058, 1079, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 698,
                                                                       713, 1079, 1100, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 713,
                                                                       728, 1100, 1121, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 728,
                                                                       743, 1121, 1142, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 743,
                                                                       758, 1142, 1163, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 758,
                                                                       773, 1163, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 773,
                                                                       788, 1184, 1205, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 788,
                                                                       803, 1205, 1226, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 803,
                                                                       818, 1226, 1247, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 848,
                                                                       869, 1268, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1808, 0, 3, 869,
                                                                       890, 1296, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1844, 0, 3, 890,
                                                                       911, 1324, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1880, 0, 3, 911,
                                                                       932, 1352, 1380, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1916, 0, 3, 932,
                                                                       953, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1952, 0, 3, 953,
                                                                       974, 1408, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 974,
                                                                       995, 1436, 1464, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 995,
                                                                       1016, 1464, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2060, 0, 3, 1058,
                                                                       1079, 1520, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2096, 0, 3, 1079,
                                                                       1100, 1548, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2132, 0, 3, 1100,
                                                                       1121, 1576, 1604, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2168, 0, 3, 1121,
                                                                       1142, 1604, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2204, 0, 3, 1142,
                                                                       1163, 1632, 1660, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2240, 0, 3, 1163,
                                                                       1184, 1660, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 1184,
                                                                       1205, 1688, 1716, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2312, 0, 3, 1205,
                                                                       1226, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2348, 0, 3, 1268,
                                                                       1296, 1772, 1808, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2393, 0, 3, 1296,
                                                                       1324, 1808, 1844, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2438, 0, 3, 1324,
                                                                       1352, 1844, 1880, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2483, 0, 3, 1352,
                                                                       1380, 1880, 1916, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2528, 0, 3, 1380,
                                                                       1408, 1916, 1952, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2573, 0, 3, 1408,
                                                                       1436, 1952, 1988, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2618, 0, 3, 1436,
                                                                       1464, 1988, 2024, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2663, 0, 3, 1520,
                                                                       1548, 2060, 2096, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2708, 0, 3, 1548,
                                                                       1576, 2096, 2132, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2753, 0, 3, 1576,
                                                                       1604, 2132, 2168, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2798, 0, 3, 1604,
                                                                       1632, 2168, 2204, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2843, 0, 3, 1632,
                                                                       1660, 2204, 2240, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2888, 0, 3, 1660,
                                                                       1688, 2240, 2276, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 2933, 0, 3, 1688,
                                                                       1716, 2276, 2312, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2978, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2981, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2984, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2987, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2990, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2993, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2996, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2999, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3002, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3005, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3008, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3011, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3014, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3017, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3020, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3023, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3026, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3029, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3032, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3035, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3038, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3041, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3044, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3047, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3050, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3053, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3056, 3, 9, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3065, 3, 10, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3074, 3, 11, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3083, 3, 12, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3092, 3, 13, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3101, 3, 14, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3110, 3, 15, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3119, 3, 16, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3128, 3, 17, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3137, 3, 18, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3146, 3, 19, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3155, 3, 20, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3164, 3, 25, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3173, 3, 26, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3182, 3, 27, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3191, 3, 28, 95,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3200, 3, 29, 98,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3209, 3, 30, 101,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3218, 3, 31, 104,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3227, 3, 32, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3236, 3, 33, 110,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3245, 3, 34, 113,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3254, 3, 35, 116,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3263, 3, 36, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3272, 3, 44, 134,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3290, 3, 47, 140,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3308, 3, 50, 146,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3326, 3, 53, 152,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3344, 3, 56, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3362, 3, 59, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3380, 3, 62, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3398, 3, 65, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3416, 3, 68, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3434, 3, 71, 188,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3452, 3, 74, 194,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3470, 3, 86, 212,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3488, 3, 89, 218,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3506, 3, 92, 224,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3524, 3, 95, 230,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3542, 3, 98, 236,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3560, 3, 101, 242,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3578, 3, 104, 248,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3596, 3, 107, 254,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3614, 3, 110, 260,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3632, 3, 113, 266,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3650, 3, 116, 272,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3668, 3, 134, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3698, 3, 140, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3728, 3, 146, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3758, 3, 152, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3788, 3, 158, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3818, 3, 164, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3848, 3, 170, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3878, 3, 176, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3908, 3, 182, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3938, 3, 188, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3968, 3, 212, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3998, 3, 218, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4028, 3, 224, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4058, 3, 230, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4088, 3, 236, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4118, 3, 242, 468,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4148, 3, 248, 478,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4178, 3, 254, 488,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4208, 3, 260, 498,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4238, 3, 266, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4268, 3, 298, 548,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4313, 3, 308, 563,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4358, 3, 318, 578,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4403, 3, 328, 593,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4448, 3, 338, 608,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4493, 3, 348, 623,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4538, 3, 358, 638,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4583, 3, 368, 653,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4628, 3, 378, 668,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4673, 3, 418, 713,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4718, 3, 428, 728,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4763, 3, 438, 743,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4808, 3, 448, 758,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4853, 3, 458, 773,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4898, 3, 468, 788,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4943, 3, 478, 803,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4988, 3, 488, 818,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5033, 3, 498, 833,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5078, 3, 548, 890,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5141, 3, 563, 911,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5204, 3, 578, 932,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5267, 3, 593, 953,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5330, 3, 608, 974,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5393, 3, 623, 995,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5456, 3, 638,
                                                                       1016, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5519, 3, 653,
                                                                       1037, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5582, 3, 713,
                                                                       1100, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5645, 3, 728,
                                                                       1121, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5708, 3, 743,
                                                                       1142, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5771, 3, 758,
                                                                       1163, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5834, 3, 773,
                                                                       1184, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5897, 3, 788,
                                                                       1205, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 5960, 3, 803,
                                                                       1226, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 6023, 3, 818,
                                                                       1247, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6086, 3, 890,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6170, 3, 911,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6254, 3, 932,
                                                                       1380, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6338, 3, 953,
                                                                       1408, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6422, 3, 974,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6506, 3, 995,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6590, 3, 1016,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6674, 3, 1100,
                                                                       1576, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6758, 3, 1121,
                                                                       1604, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6842, 3, 1142,
                                                                       1632, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 6926, 3, 1163,
                                                                       1660, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7010, 3, 1184,
                                                                       1688, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7094, 3, 1205,
                                                                       1716, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 7178, 3, 1226,
                                                                       1744, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7262, 3, 1324,
                                                                       1844, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7370, 3, 1352,
                                                                       1880, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7478, 3, 1380,
                                                                       1916, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7586, 3, 1408,
                                                                       1952, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7694, 3, 1436,
                                                                       1988, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7802, 3, 1464,
                                                                       2024, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 7910, 3, 1576,
                                                                       2132, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8018, 3, 1604,
                                                                       2168, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8126, 3, 1632,
                                                                       2204, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8234, 3, 1660,
                                                                       2240, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8342, 3, 1688,
                                                                       2276, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 8450, 3, 1716,
                                                                       2312, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8558, 3, 1844,
                                                                       2438, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8693, 3, 1880,
                                                                       2483, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8828, 3, 1916,
                                                                       2528, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 8963, 3, 1952,
                                                                       2573, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9098, 3, 1988,
                                                                       2618, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9233, 3, 2132,
                                                                       2753, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9368, 3, 2168,
                                                                       2798, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9503, 3, 2204,
                                                                       2843, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9638, 3, 2240,
                                                                       2888, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 9773, 3, 2276,
                                                                       2933, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9908, 3, 7, 8,
                                                                       2978, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9914, 3, 8, 9,
                                                                       2981, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9920, 3, 9, 10,
                                                                       2984, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9926, 3, 10, 11,
                                                                       2987, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9932, 3, 11, 12,
                                                                       2990, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9938, 3, 12, 13,
                                                                       2993, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9944, 3, 13, 14,
                                                                       2996, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9950, 3, 14, 15,
                                                                       2999, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9956, 3, 15, 16,
                                                                       3002, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9962, 3, 16, 17,
                                                                       3005, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9968, 3, 17, 18,
                                                                       3008, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9974, 3, 18, 19,
                                                                       3011, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9980, 3, 19, 20,
                                                                       3014, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9986, 3, 23, 24,
                                                                       3017, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9992, 3, 24, 25,
                                                                       3020, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9998, 3, 25, 26,
                                                                       3023, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10004, 3, 26, 27,
                                                                       3026, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10010, 3, 27, 28,
                                                                       3029, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10016, 3, 28, 29,
                                                                       3032, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10022, 3, 29, 30,
                                                                       3035, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10028, 3, 30, 31,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10034, 3, 31, 32,
                                                                       3041, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10040, 3, 32, 33,
                                                                       3044, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10046, 3, 33, 34,
                                                                       3047, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10052, 3, 34, 35,
                                                                       3050, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10058, 3, 35, 36,
                                                                       3053, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10064, 0, 3, 9908,
                                                                       2978, 9914, 3056, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10082, 0, 3, 9914,
                                                                       2981, 9920, 3065, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10100, 0, 3, 9920,
                                                                       2984, 9926, 3074, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10118, 0, 3, 9926,
                                                                       2987, 9932, 3083, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10136, 0, 3, 9932,
                                                                       2990, 9938, 3092, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10154, 0, 3, 9938,
                                                                       2993, 9944, 3101, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10172, 0, 3, 9944,
                                                                       2996, 9950, 3110, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10190, 0, 3, 9950,
                                                                       2999, 9956, 3119, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10208, 0, 3, 9956,
                                                                       3002, 9962, 3128, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10226, 0, 3, 9962,
                                                                       3005, 9968, 3137, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10244, 0, 3, 9968,
                                                                       3008, 9974, 3146, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10262, 0, 3, 9974,
                                                                       3011, 9980, 3155, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10280, 0, 3, 9986,
                                                                       3017, 9992, 3164, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10298, 0, 3, 9992,
                                                                       3020, 9998, 3173, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10316, 0, 3, 9998,
                                                                       3023, 10004, 3182, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10334, 0, 3,
                                                                       10004, 3026, 10010, 3191,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10352, 0, 3,
                                                                       10010, 3029, 10016, 3200,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10370, 0, 3,
                                                                       10016, 3032, 10022, 3209,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10388, 0, 3,
                                                                       10022, 3035, 10028, 3218,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10406, 0, 3,
                                                                       10028, 3038, 10034, 3227,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10424, 0, 3,
                                                                       10034, 3041, 10040, 3236,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10442, 0, 3,
                                                                       10040, 3044, 10046, 3245,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10460, 0, 3,
                                                                       10046, 3047, 10052, 3254,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10478, 0, 3,
                                                                       10052, 3050, 10058, 3263,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10496, 0, 3,
                                                                       10064, 3056, 10082, 122,
                                                                       128, 3272, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10532, 0, 3,
                                                                       10082, 3065, 10100, 128,
                                                                       134, 3290, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10568, 0, 3,
                                                                       10100, 3074, 10118, 134,
                                                                       140, 3308, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10604, 0, 3,
                                                                       10118, 3083, 10136, 140,
                                                                       146, 3326, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10640, 0, 3,
                                                                       10136, 3092, 10154, 146,
                                                                       152, 3344, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10676, 0, 3,
                                                                       10154, 3101, 10172, 152,
                                                                       158, 3362, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10712, 0, 3,
                                                                       10172, 3110, 10190, 158,
                                                                       164, 3380, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10748, 0, 3,
                                                                       10190, 3119, 10208, 164,
                                                                       170, 3398, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10784, 0, 3,
                                                                       10208, 3128, 10226, 170,
                                                                       176, 3416, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10820, 0, 3,
                                                                       10226, 3137, 10244, 176,
                                                                       182, 3434, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10856, 0, 3,
                                                                       10244, 3146, 10262, 182,
                                                                       188, 3452, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10892, 0, 3,
                                                                       10280, 3164, 10298, 200,
                                                                       206, 3470, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10928, 0, 3,
                                                                       10298, 3173, 10316, 206,
                                                                       212, 3488, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 10964, 0, 3,
                                                                       10316, 3182, 10334, 212,
                                                                       218, 3506, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11000, 0, 3,
                                                                       10334, 3191, 10352, 218,
                                                                       224, 3524, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11036, 0, 3,
                                                                       10352, 3200, 10370, 224,
                                                                       230, 3542, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11072, 0, 3,
                                                                       10370, 3209, 10388, 230,
                                                                       236, 3560, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11108, 0, 3,
                                                                       10388, 3218, 10406, 236,
                                                                       242, 3578, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11144, 0, 3,
                                                                       10406, 3227, 10424, 242,
                                                                       248, 3596, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11180, 0, 3,
                                                                       10424, 3236, 10442, 248,
                                                                       254, 3614, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11216, 0, 3,
                                                                       10442, 3245, 10460, 254,
                                                                       260, 3632, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11252, 0, 3,
                                                                       10460, 3254, 10478, 260,
                                                                       266, 3650, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11288, 0, 3,
                                                                       10496, 3272, 10532, 278,
                                                                       288, 3668, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11348, 0, 3,
                                                                       10532, 3290, 10568, 288,
                                                                       298, 3698, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11408, 0, 3,
                                                                       10568, 3308, 10604, 298,
                                                                       308, 3728, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11468, 0, 3,
                                                                       10604, 3326, 10640, 308,
                                                                       318, 3758, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11528, 0, 3,
                                                                       10640, 3344, 10676, 318,
                                                                       328, 3788, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11588, 0, 3,
                                                                       10676, 3362, 10712, 328,
                                                                       338, 3818, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11648, 0, 3,
                                                                       10712, 3380, 10748, 338,
                                                                       348, 3848, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11708, 0, 3,
                                                                       10748, 3398, 10784, 348,
                                                                       358, 3878, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11768, 0, 3,
                                                                       10784, 3416, 10820, 358,
                                                                       368, 3908, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11828, 0, 3,
                                                                       10820, 3434, 10856, 368,
                                                                       378, 3938, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11888, 0, 3,
                                                                       10892, 3470, 10928, 398,
                                                                       408, 3968, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11948, 0, 3,
                                                                       10928, 3488, 10964, 408,
                                                                       418, 3998, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12008, 0, 3,
                                                                       10964, 3506, 11000, 418,
                                                                       428, 4028, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12068, 0, 3,
                                                                       11000, 3524, 11036, 428,
                                                                       438, 4058, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12128, 0, 3,
                                                                       11036, 3542, 11072, 438,
                                                                       448, 4088, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12188, 0, 3,
                                                                       11072, 3560, 11108, 448,
                                                                       458, 4118, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12248, 0, 3,
                                                                       11108, 3578, 11144, 458,
                                                                       468, 4148, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12308, 0, 3,
                                                                       11144, 3596, 11180, 468,
                                                                       478, 4178, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12368, 0, 3,
                                                                       11180, 3614, 11216, 478,
                                                                       488, 4208, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12428, 0, 3,
                                                                       11216, 3632, 11252, 488,
                                                                       498, 4238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12488, 0, 3,
                                                                       11288, 3668, 11348, 518,
                                                                       533, 4268, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12578, 0, 3,
                                                                       11348, 3698, 11408, 533,
                                                                       548, 4313, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12668, 0, 3,
                                                                       11408, 3728, 11468, 548,
                                                                       563, 4358, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12758, 0, 3,
                                                                       11468, 3758, 11528, 563,
                                                                       578, 4403, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12848, 0, 3,
                                                                       11528, 3788, 11588, 578,
                                                                       593, 4448, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12938, 0, 3,
                                                                       11588, 3818, 11648, 593,
                                                                       608, 4493, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13028, 0, 3,
                                                                       11648, 3848, 11708, 608,
                                                                       623, 4538, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13118, 0, 3,
                                                                       11708, 3878, 11768, 623,
                                                                       638, 4583, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13208, 0, 3,
                                                                       11768, 3908, 11828, 638,
                                                                       653, 4628, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13298, 0, 3,
                                                                       11888, 3968, 11948, 683,
                                                                       698, 4673, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13388, 0, 3,
                                                                       11948, 3998, 12008, 698,
                                                                       713, 4718, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13478, 0, 3,
                                                                       12008, 4028, 12068, 713,
                                                                       728, 4763, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13568, 0, 3,
                                                                       12068, 4058, 12128, 728,
                                                                       743, 4808, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13658, 0, 3,
                                                                       12128, 4088, 12188, 743,
                                                                       758, 4853, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13748, 0, 3,
                                                                       12188, 4118, 12248, 758,
                                                                       773, 4898, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13838, 0, 3,
                                                                       12248, 4148, 12308, 773,
                                                                       788, 4943, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 13928, 0, 3,
                                                                       12308, 4178, 12368, 788,
                                                                       803, 4988, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14018, 0, 3,
                                                                       12368, 4208, 12428, 803,
                                                                       818, 5033, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14108, 0, 3,
                                                                       12488, 4268, 12578, 848,
                                                                       869, 5078, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14234, 0, 3,
                                                                       12578, 4313, 12668, 869,
                                                                       890, 5141, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14360, 0, 3,
                                                                       12668, 4358, 12758, 890,
                                                                       911, 5204, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14486, 0, 3,
                                                                       12758, 4403, 12848, 911,
                                                                       932, 5267, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14612, 0, 3,
                                                                       12848, 4448, 12938, 932,
                                                                       953, 5330, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14738, 0, 3,
                                                                       12938, 4493, 13028, 953,
                                                                       974, 5393, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14864, 0, 3,
                                                                       13028, 4538, 13118, 974,
                                                                       995, 5456, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 14990, 0, 3,
                                                                       13118, 4583, 13208, 995,
                                                                       1016, 5519, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15116, 0, 3,
                                                                       13298, 4673, 13388, 1058,
                                                                       1079, 5582, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15242, 0, 3,
                                                                       13388, 4718, 13478, 1079,
                                                                       1100, 5645, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15368, 0, 3,
                                                                       13478, 4763, 13568, 1100,
                                                                       1121, 5708, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15494, 0, 3,
                                                                       13568, 4808, 13658, 1121,
                                                                       1142, 5771, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15620, 0, 3,
                                                                       13658, 4853, 13748, 1142,
                                                                       1163, 5834, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15746, 0, 3,
                                                                       13748, 4898, 13838, 1163,
                                                                       1184, 5897, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15872, 0, 3,
                                                                       13838, 4943, 13928, 1184,
                                                                       1205, 5960, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 15998, 0, 3,
                                                                       13928, 4988, 14018, 1205,
                                                                       1226, 6023, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16124, 0, 3,
                                                                       14108, 5078, 14234, 1268,
                                                                       1296, 6086, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16292, 0, 3,
                                                                       14234, 5141, 14360, 1296,
                                                                       1324, 6170, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16460, 0, 3,
                                                                       14360, 5204, 14486, 1324,
                                                                       1352, 6254, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16628, 0, 3,
                                                                       14486, 5267, 14612, 1352,
                                                                       1380, 6338, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16796, 0, 3,
                                                                       14612, 5330, 14738, 1380,
                                                                       1408, 6422, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 16964, 0, 3,
                                                                       14738, 5393, 14864, 1408,
                                                                       1436, 6506, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17132, 0, 3,
                                                                       14864, 5456, 14990, 1436,
                                                                       1464, 6590, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17300, 0, 3,
                                                                       15116, 5582, 15242, 1520,
                                                                       1548, 6674, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17468, 0, 3,
                                                                       15242, 5645, 15368, 1548,
                                                                       1576, 6758, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17636, 0, 3,
                                                                       15368, 5708, 15494, 1576,
                                                                       1604, 6842, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17804, 0, 3,
                                                                       15494, 5771, 15620, 1604,
                                                                       1632, 6926, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 17972, 0, 3,
                                                                       15620, 5834, 15746, 1632,
                                                                       1660, 7010, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18140, 0, 3,
                                                                       15746, 5897, 15872, 1660,
                                                                       1688, 7094, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 18308, 0, 3,
                                                                       15872, 5960, 15998, 1688,
                                                                       1716, 7178, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18476, 0, 3,
                                                                       16124, 6086, 16292, 1772,
                                                                       1808, 7262, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18692, 0, 3,
                                                                       16292, 6170, 16460, 1808,
                                                                       1844, 7370, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 18908, 0, 3,
                                                                       16460, 6254, 16628, 1844,
                                                                       1880, 7478, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19124, 0, 3,
                                                                       16628, 6338, 16796, 1880,
                                                                       1916, 7586, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19340, 0, 3,
                                                                       16796, 6422, 16964, 1916,
                                                                       1952, 7694, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19556, 0, 3,
                                                                       16964, 6506, 17132, 1952,
                                                                       1988, 7802, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19772, 0, 3,
                                                                       17300, 6674, 17468, 2060,
                                                                       2096, 7910, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 19988, 0, 3,
                                                                       17468, 6758, 17636, 2096,
                                                                       2132, 8018, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20204, 0, 3,
                                                                       17636, 6842, 17804, 2132,
                                                                       2168, 8126, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20420, 0, 3,
                                                                       17804, 6926, 17972, 2168,
                                                                       2204, 8234, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20636, 0, 3,
                                                                       17972, 7010, 18140, 2204,
                                                                       2240, 8342, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 20852, 0, 3,
                                                                       18140, 7094, 18308, 2240,
                                                                       2276, 8450, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21068, 0, 3,
                                                                       18476, 7262, 18692, 2348,
                                                                       2393, 8558, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21338, 0, 3,
                                                                       18692, 7370, 18908, 2393,
                                                                       2438, 8693, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21608, 0, 3,
                                                                       18908, 7478, 19124, 2438,
                                                                       2483, 8828, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 21878, 0, 3,
                                                                       19124, 7586, 19340, 2483,
                                                                       2528, 8963, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22148, 0, 3,
                                                                       19340, 7694, 19556, 2528,
                                                                       2573, 9098, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22418, 0, 3,
                                                                       19772, 7910, 19988, 2663,
                                                                       2708, 9233, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22688, 0, 3,
                                                                       19988, 8018, 20204, 2708,
                                                                       2753, 9368, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 22958, 0, 3,
                                                                       20204, 8126, 20420, 2753,
                                                                       2798, 9503, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23228, 0, 3,
                                                                       20420, 8234, 20636, 2798,
                                                                       2843, 9638, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 23498, 0, 3,
                                                                       20636, 8342, 20852, 2843,
                                                                       2888, 9773, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23768, 3, 2978,
                                                                       2981, 9920, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23778, 3, 2981,
                                                                       2984, 9926, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23788, 3, 2984,
                                                                       2987, 9932, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23798, 3, 2987,
                                                                       2990, 9938, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23808, 3, 2990,
                                                                       2993, 9944, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23818, 3, 2993,
                                                                       2996, 9950, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23828, 3, 2996,
                                                                       2999, 9956, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23838, 3, 2999,
                                                                       3002, 9962, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23848, 3, 3002,
                                                                       3005, 9968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23858, 3, 3005,
                                                                       3008, 9974, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23868, 3, 3008,
                                                                       3011, 9980, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23878, 3, 3017,
                                                                       3020, 9998, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23888, 3, 3020,
                                                                       3023, 10004, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23898, 3, 3023,
                                                                       3026, 10010, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23908, 3, 3026,
                                                                       3029, 10016, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23918, 3, 3029,
                                                                       3032, 10022, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23928, 3, 3032,
                                                                       3035, 10028, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23938, 3, 3035,
                                                                       3038, 10034, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23948, 3, 3038,
                                                                       3041, 10040, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23958, 3, 3041,
                                                                       3044, 10046, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23968, 3, 3044,
                                                                       3047, 10052, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23978, 3, 3047,
                                                                       3050, 10058, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 23988, 0, 3,
                                                                       23768, 9920, 23778, 3056,
                                                                       3065, 10100, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24018, 0, 3,
                                                                       23778, 9926, 23788, 3065,
                                                                       3074, 10118, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24048, 0, 3,
                                                                       23788, 9932, 23798, 3074,
                                                                       3083, 10136, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24078, 0, 3,
                                                                       23798, 9938, 23808, 3083,
                                                                       3092, 10154, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24108, 0, 3,
                                                                       23808, 9944, 23818, 3092,
                                                                       3101, 10172, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24138, 0, 3,
                                                                       23818, 9950, 23828, 3101,
                                                                       3110, 10190, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24168, 0, 3,
                                                                       23828, 9956, 23838, 3110,
                                                                       3119, 10208, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24198, 0, 3,
                                                                       23838, 9962, 23848, 3119,
                                                                       3128, 10226, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24228, 0, 3,
                                                                       23848, 9968, 23858, 3128,
                                                                       3137, 10244, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24258, 0, 3,
                                                                       23858, 9974, 23868, 3137,
                                                                       3146, 10262, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24288, 0, 3,
                                                                       23878, 9998, 23888, 3164,
                                                                       3173, 10316, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24318, 0, 3,
                                                                       23888, 10004, 23898, 3173,
                                                                       3182, 10334, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24348, 0, 3,
                                                                       23898, 10010, 23908, 3182,
                                                                       3191, 10352, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24378, 0, 3,
                                                                       23908, 10016, 23918, 3191,
                                                                       3200, 10370, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24408, 0, 3,
                                                                       23918, 10022, 23928, 3200,
                                                                       3209, 10388, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24438, 0, 3,
                                                                       23928, 10028, 23938, 3209,
                                                                       3218, 10406, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24468, 0, 3,
                                                                       23938, 10034, 23948, 3218,
                                                                       3227, 10424, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24498, 0, 3,
                                                                       23948, 10040, 23958, 3227,
                                                                       3236, 10442, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24528, 0, 3,
                                                                       23958, 10046, 23968, 3236,
                                                                       3245, 10460, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 24558, 0, 3,
                                                                       23968, 10052, 23978, 3245,
                                                                       3254, 10478, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24588, 0, 3,
                                                                       23988, 10100, 24018, 3272,
                                                                       3290, 10568, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24648, 0, 3,
                                                                       24018, 10118, 24048, 3290,
                                                                       3308, 10604, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24708, 0, 3,
                                                                       24048, 10136, 24078, 3308,
                                                                       3326, 10640, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24768, 0, 3,
                                                                       24078, 10154, 24108, 3326,
                                                                       3344, 10676, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24828, 0, 3,
                                                                       24108, 10172, 24138, 3344,
                                                                       3362, 10712, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24888, 0, 3,
                                                                       24138, 10190, 24168, 3362,
                                                                       3380, 10748, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 24948, 0, 3,
                                                                       24168, 10208, 24198, 3380,
                                                                       3398, 10784, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25008, 0, 3,
                                                                       24198, 10226, 24228, 3398,
                                                                       3416, 10820, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25068, 0, 3,
                                                                       24228, 10244, 24258, 3416,
                                                                       3434, 10856, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25128, 0, 3,
                                                                       24288, 10316, 24318, 3470,
                                                                       3488, 10964, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25188, 0, 3,
                                                                       24318, 10334, 24348, 3488,
                                                                       3506, 11000, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25248, 0, 3,
                                                                       24348, 10352, 24378, 3506,
                                                                       3524, 11036, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25308, 0, 3,
                                                                       24378, 10370, 24408, 3524,
                                                                       3542, 11072, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25368, 0, 3,
                                                                       24408, 10388, 24438, 3542,
                                                                       3560, 11108, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25428, 0, 3,
                                                                       24438, 10406, 24468, 3560,
                                                                       3578, 11144, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25488, 0, 3,
                                                                       24468, 10424, 24498, 3578,
                                                                       3596, 11180, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25548, 0, 3,
                                                                       24498, 10442, 24528, 3596,
                                                                       3614, 11216, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 25608, 0, 3,
                                                                       24528, 10460, 24558, 3614,
                                                                       3632, 11252, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25668, 0, 3,
                                                                       24588, 10568, 24648, 3668,
                                                                       3698, 11408, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25768, 0, 3,
                                                                       24648, 10604, 24708, 3698,
                                                                       3728, 11468, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25868, 0, 3,
                                                                       24708, 10640, 24768, 3728,
                                                                       3758, 11528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 25968, 0, 3,
                                                                       24768, 10676, 24828, 3758,
                                                                       3788, 11588, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26068, 0, 3,
                                                                       24828, 10712, 24888, 3788,
                                                                       3818, 11648, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26168, 0, 3,
                                                                       24888, 10748, 24948, 3818,
                                                                       3848, 11708, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26268, 0, 3,
                                                                       24948, 10784, 25008, 3848,
                                                                       3878, 11768, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26368, 0, 3,
                                                                       25008, 10820, 25068, 3878,
                                                                       3908, 11828, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26468, 0, 3,
                                                                       25128, 10964, 25188, 3968,
                                                                       3998, 12008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26568, 0, 3,
                                                                       25188, 11000, 25248, 3998,
                                                                       4028, 12068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26668, 0, 3,
                                                                       25248, 11036, 25308, 4028,
                                                                       4058, 12128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26768, 0, 3,
                                                                       25308, 11072, 25368, 4058,
                                                                       4088, 12188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26868, 0, 3,
                                                                       25368, 11108, 25428, 4088,
                                                                       4118, 12248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 26968, 0, 3,
                                                                       25428, 11144, 25488, 4118,
                                                                       4148, 12308, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 27068, 0, 3,
                                                                       25488, 11180, 25548, 4148,
                                                                       4178, 12368, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 27168, 0, 3,
                                                                       25548, 11216, 25608, 4178,
                                                                       4208, 12428, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27268, 0, 3,
                                                                       25668, 11408, 25768, 4268,
                                                                       4313, 12668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27418, 0, 3,
                                                                       25768, 11468, 25868, 4313,
                                                                       4358, 12758, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27568, 0, 3,
                                                                       25868, 11528, 25968, 4358,
                                                                       4403, 12848, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27718, 0, 3,
                                                                       25968, 11588, 26068, 4403,
                                                                       4448, 12938, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27868, 0, 3,
                                                                       26068, 11648, 26168, 4448,
                                                                       4493, 13028, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 28018, 0, 3,
                                                                       26168, 11708, 26268, 4493,
                                                                       4538, 13118, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 28168, 0, 3,
                                                                       26268, 11768, 26368, 4538,
                                                                       4583, 13208, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 28318, 0, 3,
                                                                       26468, 12008, 26568, 4673,
                                                                       4718, 13478, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 28468, 0, 3,
                                                                       26568, 12068, 26668, 4718,
                                                                       4763, 13568, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 28618, 0, 3,
                                                                       26668, 12128, 26768, 4763,
                                                                       4808, 13658, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 28768, 0, 3,
                                                                       26768, 12188, 26868, 4808,
                                                                       4853, 13748, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 28918, 0, 3,
                                                                       26868, 12248, 26968, 4853,
                                                                       4898, 13838, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 29068, 0, 3,
                                                                       26968, 12308, 27068, 4898,
                                                                       4943, 13928, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 29218, 0, 3,
                                                                       27068, 12368, 27168, 4943,
                                                                       4988, 14018, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29368, 0, 3,
                                                                       27268, 12668, 27418, 5078,
                                                                       5141, 14360, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29578, 0, 3,
                                                                       27418, 12758, 27568, 5141,
                                                                       5204, 14486, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29788, 0, 3,
                                                                       27568, 12848, 27718, 5204,
                                                                       5267, 14612, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 29998, 0, 3,
                                                                       27718, 12938, 27868, 5267,
                                                                       5330, 14738, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 30208, 0, 3,
                                                                       27868, 13028, 28018, 5330,
                                                                       5393, 14864, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 30418, 0, 3,
                                                                       28018, 13118, 28168, 5393,
                                                                       5456, 14990, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 30628, 0, 3,
                                                                       28318, 13478, 28468, 5582,
                                                                       5645, 15368, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 30838, 0, 3,
                                                                       28468, 13568, 28618, 5645,
                                                                       5708, 15494, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 31048, 0, 3,
                                                                       28618, 13658, 28768, 5708,
                                                                       5771, 15620, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 31258, 0, 3,
                                                                       28768, 13748, 28918, 5771,
                                                                       5834, 15746, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 31468, 0, 3,
                                                                       28918, 13838, 29068, 5834,
                                                                       5897, 15872, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 31678, 0, 3,
                                                                       29068, 13928, 29218, 5897,
                                                                       5960, 15998, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 31888, 0, 3,
                                                                       29368, 14360, 29578, 6086,
                                                                       6170, 16460, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 32168, 0, 3,
                                                                       29578, 14486, 29788, 6170,
                                                                       6254, 16628, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 32448, 0, 3,
                                                                       29788, 14612, 29998, 6254,
                                                                       6338, 16796, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 32728, 0, 3,
                                                                       29998, 14738, 30208, 6338,
                                                                       6422, 16964, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 33008, 0, 3,
                                                                       30208, 14864, 30418, 6422,
                                                                       6506, 17132, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 33288, 0, 3,
                                                                       30628, 15368, 30838, 6674,
                                                                       6758, 17636, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 33568, 0, 3,
                                                                       30838, 15494, 31048, 6758,
                                                                       6842, 17804, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 33848, 0, 3,
                                                                       31048, 15620, 31258, 6842,
                                                                       6926, 17972, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34128, 0, 3,
                                                                       31258, 15746, 31468, 6926,
                                                                       7010, 18140, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 34408, 0, 3,
                                                                       31468, 15872, 31678, 7010,
                                                                       7094, 18308, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 34688, 0, 3,
                                                                       31888, 16460, 32168, 7262,
                                                                       7370, 18908, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 35048, 0, 3,
                                                                       32168, 16628, 32448, 7370,
                                                                       7478, 19124, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 35408, 0, 3,
                                                                       32448, 16796, 32728, 7478,
                                                                       7586, 19340, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 35768, 0, 3,
                                                                       32728, 16964, 33008, 7586,
                                                                       7694, 19556, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36128, 0, 3,
                                                                       33288, 17636, 33568, 7910,
                                                                       8018, 20204, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36488, 0, 3,
                                                                       33568, 17804, 33848, 8018,
                                                                       8126, 20420, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 36848, 0, 3,
                                                                       33848, 17972, 34128, 8126,
                                                                       8234, 20636, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 37208, 0, 3,
                                                                       34128, 18140, 34408, 8234,
                                                                       8342, 20852, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 37568, 0, 3,
                                                                       34688, 18908, 35048, 8558,
                                                                       8693, 21608, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 38018, 0, 3,
                                                                       35048, 19124, 35408, 8693,
                                                                       8828, 21878, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 38468, 0, 3,
                                                                       35408, 19340, 35768, 8828,
                                                                       8963, 22148, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 38918, 0, 3,
                                                                       36128, 20204, 36488, 9233,
                                                                       9368, 22958, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 39368, 0, 3,
                                                                       36488, 20420, 36848, 9368,
                                                                       9503, 23228, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 39818, 0, 3,
                                                                       36848, 20636, 37208, 9503,
                                                                       9638, 23498, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40268, 3, 9908,
                                                                       9914, 23768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40283, 3, 9914,
                                                                       9920, 23778, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40298, 3, 9920,
                                                                       9926, 23788, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40313, 3, 9926,
                                                                       9932, 23798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40328, 3, 9932,
                                                                       9938, 23808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40343, 3, 9938,
                                                                       9944, 23818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40358, 3, 9944,
                                                                       9950, 23828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40373, 3, 9950,
                                                                       9956, 23838, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40388, 3, 9956,
                                                                       9962, 23848, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40403, 3, 9962,
                                                                       9968, 23858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40418, 3, 9968,
                                                                       9974, 23868, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40433, 3, 9986,
                                                                       9992, 23878, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40448, 3, 9992,
                                                                       9998, 23888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40463, 3, 9998,
                                                                       10004, 23898, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40478, 3, 10004,
                                                                       10010, 23908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40493, 3, 10010,
                                                                       10016, 23918, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40508, 3, 10016,
                                                                       10022, 23928, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40523, 3, 10022,
                                                                       10028, 23938, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40538, 3, 10028,
                                                                       10034, 23948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40553, 3, 10034,
                                                                       10040, 23958, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40568, 3, 10040,
                                                                       10046, 23968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 40583, 3, 10046,
                                                                       10052, 23978, ncols,
                                                                       gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 40598, 0, 3,
                                                                       40268, 23768, 40283,
                                                                       10064, 10082, 23988,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 40643, 0, 3,
                                                                       40283, 23778, 40298,
                                                                       10082, 10100, 24018,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 40688, 0, 3,
                                                                       40298, 23788, 40313,
                                                                       10100, 10118, 24048,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 40733, 0, 3,
                                                                       40313, 23798, 40328,
                                                                       10118, 10136, 24078,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 40778, 0, 3,
                                                                       40328, 23808, 40343,
                                                                       10136, 10154, 24108,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 40823, 0, 3,
                                                                       40343, 23818, 40358,
                                                                       10154, 10172, 24138,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 40868, 0, 3,
                                                                       40358, 23828, 40373,
                                                                       10172, 10190, 24168,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 40913, 0, 3,
                                                                       40373, 23838, 40388,
                                                                       10190, 10208, 24198,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 40958, 0, 3,
                                                                       40388, 23848, 40403,
                                                                       10208, 10226, 24228,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 41003, 0, 3,
                                                                       40403, 23858, 40418,
                                                                       10226, 10244, 24258,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 41048, 0, 3,
                                                                       40433, 23878, 40448,
                                                                       10280, 10298, 24288,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 41093, 0, 3,
                                                                       40448, 23888, 40463,
                                                                       10298, 10316, 24318,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 41138, 0, 3,
                                                                       40463, 23898, 40478,
                                                                       10316, 10334, 24348,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 41183, 0, 3,
                                                                       40478, 23908, 40493,
                                                                       10334, 10352, 24378,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 41228, 0, 3,
                                                                       40493, 23918, 40508,
                                                                       10352, 10370, 24408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 41273, 0, 3,
                                                                       40508, 23928, 40523,
                                                                       10370, 10388, 24438,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 41318, 0, 3,
                                                                       40523, 23938, 40538,
                                                                       10388, 10406, 24468,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 41363, 0, 3,
                                                                       40538, 23948, 40553,
                                                                       10406, 10424, 24498,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 41408, 0, 3,
                                                                       40553, 23958, 40568,
                                                                       10424, 10442, 24528,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 41453, 0, 3,
                                                                       40568, 23968, 40583,
                                                                       10442, 10460, 24558,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 41498, 0, 3,
                                                                       40598, 23988, 40643,
                                                                       10496, 10532, 24588,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 41588, 0, 3,
                                                                       40643, 24018, 40688,
                                                                       10532, 10568, 24648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 41678, 0, 3,
                                                                       40688, 24048, 40733,
                                                                       10568, 10604, 24708,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 41768, 0, 3,
                                                                       40733, 24078, 40778,
                                                                       10604, 10640, 24768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 41858, 0, 3,
                                                                       40778, 24108, 40823,
                                                                       10640, 10676, 24828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 41948, 0, 3,
                                                                       40823, 24138, 40868,
                                                                       10676, 10712, 24888,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 42038, 0, 3,
                                                                       40868, 24168, 40913,
                                                                       10712, 10748, 24948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 42128, 0, 3,
                                                                       40913, 24198, 40958,
                                                                       10748, 10784, 25008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 42218, 0, 3,
                                                                       40958, 24228, 41003,
                                                                       10784, 10820, 25068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 42308, 0, 3,
                                                                       41048, 24288, 41093,
                                                                       10892, 10928, 25128,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 42398, 0, 3,
                                                                       41093, 24318, 41138,
                                                                       10928, 10964, 25188,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 42488, 0, 3,
                                                                       41138, 24348, 41183,
                                                                       10964, 11000, 25248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 42578, 0, 3,
                                                                       41183, 24378, 41228,
                                                                       11000, 11036, 25308,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 42668, 0, 3,
                                                                       41228, 24408, 41273,
                                                                       11036, 11072, 25368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 42758, 0, 3,
                                                                       41273, 24438, 41318,
                                                                       11072, 11108, 25428,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 42848, 0, 3,
                                                                       41318, 24468, 41363,
                                                                       11108, 11144, 25488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 42938, 0, 3,
                                                                       41363, 24498, 41408,
                                                                       11144, 11180, 25548,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 43028, 0, 3,
                                                                       41408, 24528, 41453,
                                                                       11180, 11216, 25608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 43118, 0, 3,
                                                                       41498, 24588, 41588,
                                                                       11288, 11348, 25668,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 43268, 0, 3,
                                                                       41588, 24648, 41678,
                                                                       11348, 11408, 25768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 43418, 0, 3,
                                                                       41678, 24708, 41768,
                                                                       11408, 11468, 25868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 43568, 0, 3,
                                                                       41768, 24768, 41858,
                                                                       11468, 11528, 25968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 43718, 0, 3,
                                                                       41858, 24828, 41948,
                                                                       11528, 11588, 26068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 43868, 0, 3,
                                                                       41948, 24888, 42038,
                                                                       11588, 11648, 26168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 44018, 0, 3,
                                                                       42038, 24948, 42128,
                                                                       11648, 11708, 26268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 44168, 0, 3,
                                                                       42128, 25008, 42218,
                                                                       11708, 11768, 26368,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 44318, 0, 3,
                                                                       42308, 25128, 42398,
                                                                       11888, 11948, 26468,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 44468, 0, 3,
                                                                       42398, 25188, 42488,
                                                                       11948, 12008, 26568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 44618, 0, 3,
                                                                       42488, 25248, 42578,
                                                                       12008, 12068, 26668,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 44768, 0, 3,
                                                                       42578, 25308, 42668,
                                                                       12068, 12128, 26768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 44918, 0, 3,
                                                                       42668, 25368, 42758,
                                                                       12128, 12188, 26868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 45068, 0, 3,
                                                                       42758, 25428, 42848,
                                                                       12188, 12248, 26968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 45218, 0, 3,
                                                                       42848, 25488, 42938,
                                                                       12248, 12308, 27068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 45368, 0, 3,
                                                                       42938, 25548, 43028,
                                                                       12308, 12368, 27168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 45518, 0, 3,
                                                                       43118, 25668, 43268,
                                                                       12488, 12578, 27268,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 45743, 0, 3,
                                                                       43268, 25768, 43418,
                                                                       12578, 12668, 27418,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 45968, 0, 3,
                                                                       43418, 25868, 43568,
                                                                       12668, 12758, 27568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 46193, 0, 3,
                                                                       43568, 25968, 43718,
                                                                       12758, 12848, 27718,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 46418, 0, 3,
                                                                       43718, 26068, 43868,
                                                                       12848, 12938, 27868,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 46643, 0, 3,
                                                                       43868, 26168, 44018,
                                                                       12938, 13028, 28018,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 46868, 0, 3,
                                                                       44018, 26268, 44168,
                                                                       13028, 13118, 28168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 47093, 0, 3,
                                                                       44318, 26468, 44468,
                                                                       13298, 13388, 28318,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 47318, 0, 3,
                                                                       44468, 26568, 44618,
                                                                       13388, 13478, 28468,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 47543, 0, 3,
                                                                       44618, 26668, 44768,
                                                                       13478, 13568, 28618,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 47768, 0, 3,
                                                                       44768, 26768, 44918,
                                                                       13568, 13658, 28768,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 47993, 0, 3,
                                                                       44918, 26868, 45068,
                                                                       13658, 13748, 28918,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48218, 0, 3,
                                                                       45068, 26968, 45218,
                                                                       13748, 13838, 29068,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 48443, 0, 3,
                                                                       45218, 27068, 45368,
                                                                       13838, 13928, 29218,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 48668, 0, 3,
                                                                       45518, 27268, 45743,
                                                                       14108, 14234, 29368,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 48983, 0, 3,
                                                                       45743, 27418, 45968,
                                                                       14234, 14360, 29578,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 49298, 0, 3,
                                                                       45968, 27568, 46193,
                                                                       14360, 14486, 29788,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 49613, 0, 3,
                                                                       46193, 27718, 46418,
                                                                       14486, 14612, 29998,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 49928, 0, 3,
                                                                       46418, 27868, 46643,
                                                                       14612, 14738, 30208,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 50243, 0, 3,
                                                                       46643, 28018, 46868,
                                                                       14738, 14864, 30418,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 50558, 0, 3,
                                                                       47093, 28318, 47318,
                                                                       15116, 15242, 30628,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 50873, 0, 3,
                                                                       47318, 28468, 47543,
                                                                       15242, 15368, 30838,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51188, 0, 3,
                                                                       47543, 28618, 47768,
                                                                       15368, 15494, 31048,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51503, 0, 3,
                                                                       47768, 28768, 47993,
                                                                       15494, 15620, 31258,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 51818, 0, 3,
                                                                       47993, 28918, 48218,
                                                                       15620, 15746, 31468,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 52133, 0, 3,
                                                                       48218, 29068, 48443,
                                                                       15746, 15872, 31678,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 52448, 0, 3,
                                                                       48668, 29368, 48983,
                                                                       16124, 16292, 31888,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 52868, 0, 3,
                                                                       48983, 29578, 49298,
                                                                       16292, 16460, 32168,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 53288, 0, 3,
                                                                       49298, 29788, 49613,
                                                                       16460, 16628, 32448,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 53708, 0, 3,
                                                                       49613, 29998, 49928,
                                                                       16628, 16796, 32728,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 54128, 0, 3,
                                                                       49928, 30208, 50243,
                                                                       16796, 16964, 33008,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 54548, 0, 3,
                                                                       50558, 30628, 50873,
                                                                       17300, 17468, 33288,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 54968, 0, 3,
                                                                       50873, 30838, 51188,
                                                                       17468, 17636, 33568,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 55388, 0, 3,
                                                                       51188, 31048, 51503,
                                                                       17636, 17804, 33848,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 55808, 0, 3,
                                                                       51503, 31258, 51818,
                                                                       17804, 17972, 34128,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 56228, 0, 3,
                                                                       51818, 31468, 52133,
                                                                       17972, 18140, 34408,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 56648, 0, 3,
                                                                       52448, 31888, 52868,
                                                                       18476, 18692, 34688,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 57188, 0, 3,
                                                                       52868, 32168, 53288,
                                                                       18692, 18908, 35048,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 57728, 0, 3,
                                                                       53288, 32448, 53708,
                                                                       18908, 19124, 35408,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 58268, 0, 3,
                                                                       53708, 32728, 54128,
                                                                       19124, 19340, 35768,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 58808, 0, 3,
                                                                       54548, 33288, 54968,
                                                                       19772, 19988, 36128,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 59348, 0, 3,
                                                                       54968, 33568, 55388,
                                                                       19988, 20204, 36488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 59888, 0, 3,
                                                                       55388, 33848, 55808,
                                                                       20204, 20420, 36848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 60428, 0, 3,
                                                                       55808, 34128, 56228,
                                                                       20420, 20636, 37208,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 60968, 0, 3,
                                                                       56648, 34688, 57188,
                                                                       21068, 21338, 37568,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 61643, 0, 3,
                                                                       57188, 35048, 57728,
                                                                       21338, 21608, 38018,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 62318, 0, 3,
                                                                       57728, 35408, 58268,
                                                                       21608, 21878, 38468,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 62993, 0, 3,
                                                                       58808, 36128, 59348,
                                                                       22418, 22688, 38918,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 63668, 0, 3,
                                                                       59348, 36488, 59888,
                                                                       22688, 22958, 39368,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 64343, 0, 3,
                                                                       59888, 36848, 60428,
                                                                       22958, 23228, 39818,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65018, 3, 23768,
                                                                       23778, 40298, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65039, 3, 23778,
                                                                       23788, 40313, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65060, 3, 23788,
                                                                       23798, 40328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65081, 3, 23798,
                                                                       23808, 40343, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65102, 3, 23808,
                                                                       23818, 40358, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65123, 3, 23818,
                                                                       23828, 40373, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65144, 3, 23828,
                                                                       23838, 40388, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65165, 3, 23838,
                                                                       23848, 40403, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65186, 3, 23848,
                                                                       23858, 40418, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65207, 3, 23878,
                                                                       23888, 40463, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65228, 3, 23888,
                                                                       23898, 40478, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65249, 3, 23898,
                                                                       23908, 40493, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65270, 3, 23908,
                                                                       23918, 40508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65291, 3, 23918,
                                                                       23928, 40523, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65312, 3, 23928,
                                                                       23938, 40538, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65333, 3, 23938,
                                                                       23948, 40553, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65354, 3, 23948,
                                                                       23958, 40568, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65375, 3, 23958,
                                                                       23968, 40583, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65396, 0, 3,
                                                                       65018, 40298, 65039,
                                                                       23988, 24018, 40688,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65459, 0, 3,
                                                                       65039, 40313, 65060,
                                                                       24018, 24048, 40733,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65522, 0, 3,
                                                                       65060, 40328, 65081,
                                                                       24048, 24078, 40778,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65585, 0, 3,
                                                                       65081, 40343, 65102,
                                                                       24078, 24108, 40823,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65648, 0, 3,
                                                                       65102, 40358, 65123,
                                                                       24108, 24138, 40868,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65711, 0, 3,
                                                                       65123, 40373, 65144,
                                                                       24138, 24168, 40913,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65774, 0, 3,
                                                                       65144, 40388, 65165,
                                                                       24168, 24198, 40958,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65837, 0, 3,
                                                                       65165, 40403, 65186,
                                                                       24198, 24228, 41003,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65900, 0, 3,
                                                                       65207, 40463, 65228,
                                                                       24288, 24318, 41138,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 65963, 0, 3,
                                                                       65228, 40478, 65249,
                                                                       24318, 24348, 41183,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 66026, 0, 3,
                                                                       65249, 40493, 65270,
                                                                       24348, 24378, 41228,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 66089, 0, 3,
                                                                       65270, 40508, 65291,
                                                                       24378, 24408, 41273,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 66152, 0, 3,
                                                                       65291, 40523, 65312,
                                                                       24408, 24438, 41318,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 66215, 0, 3,
                                                                       65312, 40538, 65333,
                                                                       24438, 24468, 41363,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 66278, 0, 3,
                                                                       65333, 40553, 65354,
                                                                       24468, 24498, 41408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 66341, 0, 3,
                                                                       65354, 40568, 65375,
                                                                       24498, 24528, 41453,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66404, 0, 3,
                                                                       65396, 40688, 65459,
                                                                       24588, 24648, 41678,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66530, 0, 3,
                                                                       65459, 40733, 65522,
                                                                       24648, 24708, 41768,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66656, 0, 3,
                                                                       65522, 40778, 65585,
                                                                       24708, 24768, 41858,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66782, 0, 3,
                                                                       65585, 40823, 65648,
                                                                       24768, 24828, 41948,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 66908, 0, 3,
                                                                       65648, 40868, 65711,
                                                                       24828, 24888, 42038,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67034, 0, 3,
                                                                       65711, 40913, 65774,
                                                                       24888, 24948, 42128,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67160, 0, 3,
                                                                       65774, 40958, 65837,
                                                                       24948, 25008, 42218,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67286, 0, 3,
                                                                       65900, 41138, 65963,
                                                                       25128, 25188, 42488,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67412, 0, 3,
                                                                       65963, 41183, 66026,
                                                                       25188, 25248, 42578,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67538, 0, 3,
                                                                       66026, 41228, 66089,
                                                                       25248, 25308, 42668,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67664, 0, 3,
                                                                       66089, 41273, 66152,
                                                                       25308, 25368, 42758,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67790, 0, 3,
                                                                       66152, 41318, 66215,
                                                                       25368, 25428, 42848,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 67916, 0, 3,
                                                                       66215, 41363, 66278,
                                                                       25428, 25488, 42938,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 68042, 0, 3,
                                                                       66278, 41408, 66341,
                                                                       25488, 25548, 43028,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68168, 0, 3,
                                                                       66404, 41678, 66530,
                                                                       25668, 25768, 43418,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68378, 0, 3,
                                                                       66530, 41768, 66656,
                                                                       25768, 25868, 43568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68588, 0, 3,
                                                                       66656, 41858, 66782,
                                                                       25868, 25968, 43718,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 68798, 0, 3,
                                                                       66782, 41948, 66908,
                                                                       25968, 26068, 43868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 69008, 0, 3,
                                                                       66908, 42038, 67034,
                                                                       26068, 26168, 44018,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 69218, 0, 3,
                                                                       67034, 42128, 67160,
                                                                       26168, 26268, 44168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 69428, 0, 3,
                                                                       67286, 42488, 67412,
                                                                       26468, 26568, 44618,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 69638, 0, 3,
                                                                       67412, 42578, 67538,
                                                                       26568, 26668, 44768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 69848, 0, 3,
                                                                       67538, 42668, 67664,
                                                                       26668, 26768, 44918,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 70058, 0, 3,
                                                                       67664, 42758, 67790,
                                                                       26768, 26868, 45068,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 70268, 0, 3,
                                                                       67790, 42848, 67916,
                                                                       26868, 26968, 45218,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 70478, 0, 3,
                                                                       67916, 42938, 68042,
                                                                       26968, 27068, 45368,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 70688, 0, 3,
                                                                       68168, 43418, 68378,
                                                                       27268, 27418, 45968,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 71003, 0, 3,
                                                                       68378, 43568, 68588,
                                                                       27418, 27568, 46193,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 71318, 0, 3,
                                                                       68588, 43718, 68798,
                                                                       27568, 27718, 46418,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 71633, 0, 3,
                                                                       68798, 43868, 69008,
                                                                       27718, 27868, 46643,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 71948, 0, 3,
                                                                       69008, 44018, 69218,
                                                                       27868, 28018, 46868,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 72263, 0, 3,
                                                                       69428, 44618, 69638,
                                                                       28318, 28468, 47543,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 72578, 0, 3,
                                                                       69638, 44768, 69848,
                                                                       28468, 28618, 47768,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 72893, 0, 3,
                                                                       69848, 44918, 70058,
                                                                       28618, 28768, 47993,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 73208, 0, 3,
                                                                       70058, 45068, 70268,
                                                                       28768, 28918, 48218,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 73523, 0, 3,
                                                                       70268, 45218, 70478,
                                                                       28918, 29068, 48443,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 73838, 0, 3,
                                                                       70688, 45968, 71003,
                                                                       29368, 29578, 49298,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 74279, 0, 3,
                                                                       71003, 46193, 71318,
                                                                       29578, 29788, 49613,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 74720, 0, 3,
                                                                       71318, 46418, 71633,
                                                                       29788, 29998, 49928,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 75161, 0, 3,
                                                                       71633, 46643, 71948,
                                                                       29998, 30208, 50243,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 75602, 0, 3,
                                                                       72263, 47543, 72578,
                                                                       30628, 30838, 51188,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 76043, 0, 3,
                                                                       72578, 47768, 72893,
                                                                       30838, 31048, 51503,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 76484, 0, 3,
                                                                       72893, 47993, 73208,
                                                                       31048, 31258, 51818,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 76925, 0, 3,
                                                                       73208, 48218, 73523,
                                                                       31258, 31468, 52133,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 77366, 0, 3,
                                                                       73838, 49298, 74279,
                                                                       31888, 32168, 53288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 77954, 0, 3,
                                                                       74279, 49613, 74720,
                                                                       32168, 32448, 53708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 78542, 0, 3,
                                                                       74720, 49928, 75161,
                                                                       32448, 32728, 54128,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 79130, 0, 3,
                                                                       75602, 51188, 76043,
                                                                       33288, 33568, 55388,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 79718, 0, 3,
                                                                       76043, 51503, 76484,
                                                                       33568, 33848, 55808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 80306, 0, 3,
                                                                       76484, 51818, 76925,
                                                                       33848, 34128, 56228,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 80894, 0, 3,
                                                                       77366, 53288, 77954,
                                                                       34688, 35048, 57728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 81650, 0, 3,
                                                                       77954, 53708, 78542,
                                                                       35048, 35408, 58268,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 82406, 0, 3,
                                                                       79130, 55388, 79718,
                                                                       36128, 36488, 59888,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 83162, 0, 3,
                                                                       79718, 55808, 80306,
                                                                       36488, 36848, 60428,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 83918, 0, 3,
                                                                       80894, 57728, 81650,
                                                                       37568, 38018, 62318,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 84863, 0, 3,
                                                                       82406, 59888, 83162,
                                                                       38918, 39368, 64343,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85808, 3, 40268,
                                                                       40283, 65018, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85836, 3, 40283,
                                                                       40298, 65039, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85864, 3, 40298,
                                                                       40313, 65060, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85892, 3, 40313,
                                                                       40328, 65081, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85920, 3, 40328,
                                                                       40343, 65102, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85948, 3, 40343,
                                                                       40358, 65123, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85976, 3, 40358,
                                                                       40373, 65144, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86004, 3, 40373,
                                                                       40388, 65165, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86032, 3, 40388,
                                                                       40403, 65186, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86060, 3, 40433,
                                                                       40448, 65207, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86088, 3, 40448,
                                                                       40463, 65228, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86116, 3, 40463,
                                                                       40478, 65249, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86144, 3, 40478,
                                                                       40493, 65270, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86172, 3, 40493,
                                                                       40508, 65291, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86200, 3, 40508,
                                                                       40523, 65312, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86228, 3, 40523,
                                                                       40538, 65333, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86256, 3, 40538,
                                                                       40553, 65354, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 86284, 3, 40553,
                                                                       40568, 65375, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 86312, 0, 3,
                                                                       85808, 65018, 85836,
                                                                       40598, 40643, 65396,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 86396, 0, 3,
                                                                       85836, 65039, 85864,
                                                                       40643, 40688, 65459,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 86480, 0, 3,
                                                                       85864, 65060, 85892,
                                                                       40688, 40733, 65522,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 86564, 0, 3,
                                                                       85892, 65081, 85920,
                                                                       40733, 40778, 65585,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 86648, 0, 3,
                                                                       85920, 65102, 85948,
                                                                       40778, 40823, 65648,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 86732, 0, 3,
                                                                       85948, 65123, 85976,
                                                                       40823, 40868, 65711,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 86816, 0, 3,
                                                                       85976, 65144, 86004,
                                                                       40868, 40913, 65774,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 86900, 0, 3,
                                                                       86004, 65165, 86032,
                                                                       40913, 40958, 65837,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 86984, 0, 3,
                                                                       86060, 65207, 86088,
                                                                       41048, 41093, 65900,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 87068, 0, 3,
                                                                       86088, 65228, 86116,
                                                                       41093, 41138, 65963,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 87152, 0, 3,
                                                                       86116, 65249, 86144,
                                                                       41138, 41183, 66026,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 87236, 0, 3,
                                                                       86144, 65270, 86172,
                                                                       41183, 41228, 66089,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 87320, 0, 3,
                                                                       86172, 65291, 86200,
                                                                       41228, 41273, 66152,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 87404, 0, 3,
                                                                       86200, 65312, 86228,
                                                                       41273, 41318, 66215,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 87488, 0, 3,
                                                                       86228, 65333, 86256,
                                                                       41318, 41363, 66278,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 87572, 0, 3,
                                                                       86256, 65354, 86284,
                                                                       41363, 41408, 66341,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 87656, 0, 3,
                                                                       86312, 65396, 86396,
                                                                       41498, 41588, 66404,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 87824, 0, 3,
                                                                       86396, 65459, 86480,
                                                                       41588, 41678, 66530,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 87992, 0, 3,
                                                                       86480, 65522, 86564,
                                                                       41678, 41768, 66656,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 88160, 0, 3,
                                                                       86564, 65585, 86648,
                                                                       41768, 41858, 66782,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 88328, 0, 3,
                                                                       86648, 65648, 86732,
                                                                       41858, 41948, 66908,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 88496, 0, 3,
                                                                       86732, 65711, 86816,
                                                                       41948, 42038, 67034,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 88664, 0, 3,
                                                                       86816, 65774, 86900,
                                                                       42038, 42128, 67160,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 88832, 0, 3,
                                                                       86984, 65900, 87068,
                                                                       42308, 42398, 67286,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 89000, 0, 3,
                                                                       87068, 65963, 87152,
                                                                       42398, 42488, 67412,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 89168, 0, 3,
                                                                       87152, 66026, 87236,
                                                                       42488, 42578, 67538,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 89336, 0, 3,
                                                                       87236, 66089, 87320,
                                                                       42578, 42668, 67664,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 89504, 0, 3,
                                                                       87320, 66152, 87404,
                                                                       42668, 42758, 67790,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 89672, 0, 3,
                                                                       87404, 66215, 87488,
                                                                       42758, 42848, 67916,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 89840, 0, 3,
                                                                       87488, 66278, 87572,
                                                                       42848, 42938, 68042,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 90008, 0, 3,
                                                                       87656, 66404, 87824,
                                                                       43118, 43268, 68168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 90288, 0, 3,
                                                                       87824, 66530, 87992,
                                                                       43268, 43418, 68378,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 90568, 0, 3,
                                                                       87992, 66656, 88160,
                                                                       43418, 43568, 68588,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 90848, 0, 3,
                                                                       88160, 66782, 88328,
                                                                       43568, 43718, 68798,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 91128, 0, 3,
                                                                       88328, 66908, 88496,
                                                                       43718, 43868, 69008,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 91408, 0, 3,
                                                                       88496, 67034, 88664,
                                                                       43868, 44018, 69218,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 91688, 0, 3,
                                                                       88832, 67286, 89000,
                                                                       44318, 44468, 69428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 91968, 0, 3,
                                                                       89000, 67412, 89168,
                                                                       44468, 44618, 69638,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 92248, 0, 3,
                                                                       89168, 67538, 89336,
                                                                       44618, 44768, 69848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 92528, 0, 3,
                                                                       89336, 67664, 89504,
                                                                       44768, 44918, 70058,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 92808, 0, 3,
                                                                       89504, 67790, 89672,
                                                                       44918, 45068, 70268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 93088, 0, 3,
                                                                       89672, 67916, 89840,
                                                                       45068, 45218, 70478,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 93368, 0, 3,
                                                                       90008, 68168, 90288,
                                                                       45518, 45743, 70688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 93788, 0, 3,
                                                                       90288, 68378, 90568,
                                                                       45743, 45968, 71003,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 94208, 0, 3,
                                                                       90568, 68588, 90848,
                                                                       45968, 46193, 71318,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 94628, 0, 3,
                                                                       90848, 68798, 91128,
                                                                       46193, 46418, 71633,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 95048, 0, 3,
                                                                       91128, 69008, 91408,
                                                                       46418, 46643, 71948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 95468, 0, 3,
                                                                       91688, 69428, 91968,
                                                                       47093, 47318, 72263,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 95888, 0, 3,
                                                                       91968, 69638, 92248,
                                                                       47318, 47543, 72578,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 96308, 0, 3,
                                                                       92248, 69848, 92528,
                                                                       47543, 47768, 72893,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 96728, 0, 3,
                                                                       92528, 70058, 92808,
                                                                       47768, 47993, 73208,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 97148, 0, 3,
                                                                       92808, 70268, 93088,
                                                                       47993, 48218, 73523,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 97568, 0, 3,
                                                                       93368, 70688, 93788,
                                                                       48668, 48983, 73838,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 98156, 0, 3,
                                                                       93788, 71003, 94208,
                                                                       48983, 49298, 74279,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 98744, 0, 3,
                                                                       94208, 71318, 94628,
                                                                       49298, 49613, 74720,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 99332, 0, 3,
                                                                       94628, 71633, 95048,
                                                                       49613, 49928, 75161,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 99920, 0, 3,
                                                                       95468, 72263, 95888,
                                                                       50558, 50873, 75602,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 100508, 0, 3,
                                                                       95888, 72578, 96308,
                                                                       50873, 51188, 76043,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 101096, 0, 3,
                                                                       96308, 72893, 96728,
                                                                       51188, 51503, 76484,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 101684, 0, 3,
                                                                       96728, 73208, 97148,
                                                                       51503, 51818, 76925,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 102272, 0, 3,
                                                                       97568, 73838, 98156,
                                                                       52448, 52868, 77366,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 103056, 0, 3,
                                                                       98156, 74279, 98744,
                                                                       52868, 53288, 77954,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 103840, 0, 3,
                                                                       98744, 74720, 99332,
                                                                       53288, 53708, 78542,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 104624, 0, 3,
                                                                       99920, 75602, 100508,
                                                                       54548, 54968, 79130,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 105408, 0, 3,
                                                                       100508, 76043, 101096,
                                                                       54968, 55388, 79718,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 106192, 0, 3,
                                                                       101096, 76484, 101684,
                                                                       55388, 55808, 80306,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 106976, 0, 3,
                                                                       102272, 77366, 103056,
                                                                       56648, 57188, 80894,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 107984, 0, 3,
                                                                       103056, 77954, 103840,
                                                                       57188, 57728, 81650,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 108992, 0, 3,
                                                                       104624, 79130, 105408,
                                                                       58808, 59348, 82406,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 110000, 0, 3,
                                                                       105408, 79718, 106192,
                                                                       59348, 59888, 83162,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 111008, 0, 3,
                                                                       106976, 80894, 107984,
                                                                       60968, 61643, 83918,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 112268, 0, 3,
                                                                       108992, 82406, 110000,
                                                                       62993, 63668, 84863,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_f_x(buffer, 113528, 87656, 93368, 1, 28, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 113808, 87656, 93368, 1, 28, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 114088, 87656, 93368, 1, 28, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 114368, 88832, 95468, 1, 28, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 114648, 88832, 95468, 1, 28, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 114928, 88832, 95468, 1, 28, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 115208, 90008, 97568, 1, 28, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 115628, 90008, 97568, 1, 28, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 116048, 90008, 97568, 1, 28, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 116468, 91688, 99920, 1, 28, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 116888, 91688, 99920, 1, 28, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 117308, 91688, 99920, 1, 28, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 117728, 93368, 102272, 1, 28, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 118316, 93368, 102272, 1, 28, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 118904, 93368, 102272, 1, 28, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 119492, 95468, 104624, 1, 28, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 120080, 95468, 104624, 1, 28, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 120668, 95468, 104624, 1, 28, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 121256, 97568, 106976, 1, 28, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 122040, 97568, 106976, 1, 28, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 122824, 97568, 106976, 1, 28, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 123608, 99920, 108992, 1, 28, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 124392, 99920, 108992, 1, 28, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 125176, 99920, 108992, 1, 28, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 125960, 102272, 111008, 1, 28, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 126968, 102272, 111008, 1, 28, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 127976, 102272, 111008, 1, 28, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 128984, 104624, 112268, 1, 28, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 129992, 104624, 112268, 1, 28, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 131000, 104624, 112268, 1, 28, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 132008, 113528, 280, ncols);

                    simdfunc::contract_primitives(buffer, 132418, 113808, 280, ncols);

                    simdfunc::contract_primitives(buffer, 132828, 114088, 280, ncols);

                    simdfunc::contract_primitives(buffer, 133238, 90008, 280, ncols);

                    simdfunc::contract_primitives(buffer, 133648, 114368, 280, ncols);

                    simdfunc::contract_primitives(buffer, 134058, 114648, 280, ncols);

                    simdfunc::contract_primitives(buffer, 134468, 114928, 280, ncols);

                    simdfunc::contract_primitives(buffer, 134878, 91688, 280, ncols);

                    simdfunc::contract_primitives(buffer, 135288, 115208, 420, ncols);

                    simdfunc::contract_primitives(buffer, 135903, 115628, 420, ncols);

                    simdfunc::contract_primitives(buffer, 136518, 116048, 420, ncols);

                    simdfunc::contract_primitives(buffer, 137133, 93368, 420, ncols);

                    simdfunc::contract_primitives(buffer, 137748, 116468, 420, ncols);

                    simdfunc::contract_primitives(buffer, 138363, 116888, 420, ncols);

                    simdfunc::contract_primitives(buffer, 138978, 117308, 420, ncols);

                    simdfunc::contract_primitives(buffer, 139593, 95468, 420, ncols);

                    simdfunc::contract_primitives(buffer, 140208, 117728, 588, ncols);

                    simdfunc::contract_primitives(buffer, 141069, 118316, 588, ncols);

                    simdfunc::contract_primitives(buffer, 141930, 118904, 588, ncols);

                    simdfunc::contract_primitives(buffer, 142791, 97568, 588, ncols);

                    simdfunc::contract_primitives(buffer, 143652, 119492, 588, ncols);

                    simdfunc::contract_primitives(buffer, 144513, 120080, 588, ncols);

                    simdfunc::contract_primitives(buffer, 145374, 120668, 588, ncols);

                    simdfunc::contract_primitives(buffer, 146235, 99920, 588, ncols);

                    simdfunc::contract_primitives(buffer, 147096, 121256, 784, ncols);

                    simdfunc::contract_primitives(buffer, 148244, 122040, 784, ncols);

                    simdfunc::contract_primitives(buffer, 149392, 122824, 784, ncols);

                    simdfunc::contract_primitives(buffer, 150540, 102272, 784, ncols);

                    simdfunc::contract_primitives(buffer, 151688, 123608, 784, ncols);

                    simdfunc::contract_primitives(buffer, 152836, 124392, 784, ncols);

                    simdfunc::contract_primitives(buffer, 153984, 125176, 784, ncols);

                    simdfunc::contract_primitives(buffer, 155132, 104624, 784, ncols);

                    simdfunc::contract_primitives(buffer, 156280, 125960, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 157756, 126968, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 159232, 127976, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 160708, 128984, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 162184, 129992, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 163660, 131000, 1008, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 132288, 132008, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 132698, 132418, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 133108, 132828, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 133518, 133238, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 133928, 133648, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 134338, 134058, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 134748, 134468, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 135158, 134878, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 135708, 135288, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 136323, 135903, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 136938, 136518, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 137553, 137133, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 138168, 137748, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 138783, 138363, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 139398, 138978, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 140013, 139593, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 140796, 140208, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 141657, 141069, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 142518, 141930, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 143379, 142791, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 144240, 143652, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 145101, 144513, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 145962, 145374, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 146823, 146235, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 147880, 147096, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 149028, 148244, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 150176, 149392, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 151324, 150540, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 152472, 151688, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 153620, 152836, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 154768, 153984, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 155916, 155132, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 157288, 156280, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 158764, 157756, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 160240, 159232, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 161716, 160708, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 163192, 162184, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 164668, 163660, 36, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 165136, 132288,
                                                       133518, 135708, 13, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 165526, 132698,
                                                       133518, 136323, 13, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 165916, 133108,
                                                       133518, 136938, 13, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 166306, 133518, 137553, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 166696, 133928,
                                                       135158, 138168, 13, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 167086, 134338,
                                                       135158, 138783, 13, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 167476, 134748,
                                                       135158, 139398, 13, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 167866, 135158, 140013, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 168256, 135708,
                                                       137553, 140796, 13, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 168841, 136323,
                                                       137553, 141657, 13, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 169426, 136938,
                                                       137553, 142518, 13, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 170011, 137553, 143379, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 170596, 138168,
                                                       140013, 144240, 13, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 171181, 138783,
                                                       140013, 145101, 13, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 171766, 139398,
                                                       140013, 145962, 13, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 172351, 140013, 146823, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 172936, 140796,
                                                       143379, 147880, 13, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 173755, 141657,
                                                       143379, 149028, 13, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 174574, 142518,
                                                       143379, 150176, 13, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 175393, 143379, 151324, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 176212, 144240,
                                                       146823, 152472, 13, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 177031, 145101,
                                                       146823, 153620, 13, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 177850, 145962,
                                                       146823, 154768, 13, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 178669, 146823, 155916, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 179488, 147880,
                                                       151324, 157288, 13, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 180580, 149028,
                                                       151324, 158764, 13, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 181672, 150176,
                                                       151324, 160240, 13, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 182764, 152472,
                                                       155916, 161716, 13, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 183856, 153620,
                                                       155916, 163192, 13, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 184948, 154768,
                                                       155916, 164668, 13, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 186040, 165136,
                                                       166306, 168256, 13, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 186820, 165526,
                                                       166306, 168841, 13, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 187600, 165916,
                                                       166306, 169426, 13, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 188380, 166306, 170011, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 189160, 166696,
                                                       167866, 170596, 13, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 189940, 167086,
                                                       167866, 171181, 13, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 190720, 167476,
                                                       167866, 171766, 13, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 191500, 167866, 172351, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 192280, 168256,
                                                       170011, 172936, 13, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 193450, 168841,
                                                       170011, 173755, 13, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 194620, 169426,
                                                       170011, 174574, 13, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 195790, 170011, 175393, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 196960, 170596,
                                                       172351, 176212, 13, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 198130, 171181,
                                                       172351, 177031, 13, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 199300, 171766,
                                                       172351, 177850, 13, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 200470, 172351, 178669, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 201640, 172936,
                                                       175393, 179488, 13, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 203278, 173755,
                                                       175393, 180580, 13, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 204916, 174574,
                                                       175393, 181672, 13, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 206554, 176212,
                                                       178669, 182764, 13, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 208192, 177031,
                                                       178669, 183856, 13, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 209830, 177850,
                                                       178669, 184948, 13, nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 211468, 186040,
                                                       188380, 192280, 13, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 212768, 186820,
                                                       188380, 193450, 13, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 214068, 187600,
                                                       188380, 194620, 13, nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 215368, 188380, 195790, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 216668, 189160,
                                                       191500, 196960, 13, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 217968, 189940,
                                                       191500, 198130, 13, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 219268, 190720,
                                                       191500, 199300, 13, nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 220568, 191500, 200470, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 221868, 192280,
                                                       195790, 201640, 13, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 223818, 193450,
                                                       195790, 203278, 13, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 225768, 194620,
                                                       195790, 204916, 13, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 227718, 196960,
                                                       200470, 206554, 13, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 229668, 198130,
                                                       200470, 208192, 13, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 231618, 199300,
                                                       200470, 209830, 13, nmax);

        simdtrf::compute_hrr_geom_100x_fg_out_of_first(buffer, coordinates, 233568, 211468,
                                                       215368, 221868, 13, nmax);

        simdtrf::compute_hrr_geom_100y_fg_out_of_first(buffer, coordinates, 235518, 212768,
                                                       215368, 223818, 13, nmax);

        simdtrf::compute_hrr_geom_100z_fg_out_of_first(buffer, coordinates, 237468, 214068,
                                                       215368, 225768, 13, nmax);

        simdtrf::compute_hrr_geom_100x_fg_out_of_first(buffer, coordinates, 239418, 216668,
                                                       220568, 227718, 13, nmax);

        simdtrf::compute_hrr_geom_100y_fg_out_of_first(buffer, coordinates, 241368, 217968,
                                                       220568, 229668, 13, nmax);

        simdtrf::compute_hrr_geom_100z_fg_out_of_first(buffer, coordinates, 243318, 219268,
                                                       220568, 231618, 13, nmax);

        simdtrf::transform_g_inner(buffer, 245268, 239418, 10, 13, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 245268, 117, nmax);

        simdtrf::transform_g_inner(buffer, 245268, 241368, 10, 13, nmax);

        simdtrf::transform_f_outer(values + 819 * nvalues + n * npairs, nvalues, buffer, 245268,
                                   117, nmax);

        simdtrf::transform_g_inner(buffer, 245268, 243318, 10, 13, nmax);

        simdtrf::transform_f_outer(values + 1638 * nvalues + n * npairs, nvalues, buffer, 245268,
                                   117, nmax);

        simdtrf::transform_g_inner(buffer, 245268, 233568, 10, 13, nmax);

        simdtrf::transform_f_outer(values + 2457 * nvalues + n * npairs, nvalues, buffer, 245268,
                                   117, nmax);

        simdtrf::transform_g_inner(buffer, 245268, 235518, 10, 13, nmax);

        simdtrf::transform_f_outer(values + 3276 * nvalues + n * npairs, nvalues, buffer, 245268,
                                   117, nmax);

        simdtrf::transform_g_inner(buffer, 245268, 237468, 10, 13, nmax);

        simdtrf::transform_f_outer(values + 4095 * nvalues + n * npairs, nvalues, buffer, 245268,
                                   117, nmax);
    }

    for (size_t m = 0; m < 4914; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
