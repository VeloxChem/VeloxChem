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


#include "SimdThreeCenterElectronRepulsionRecHDI.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

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
#include "SimdTransferHD.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hdi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hdi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto dimensions = simdfunc::make_column_dimensions(
        a_function, b_function, c_function, npairs, coordinates,
        screenfunc::three_center_electron_repulsion_primitive_bound,
        threshold / static_cast<double>(nprims));

    const auto nmax = simdfunc::prepare_buffer(buffer, 47166, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 715 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 47166, 38767, 3017, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 13,
                                                             ncols, fj, mu, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 21, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 24, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 27, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 60, 0, 3, 7, 8,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 66, 0, 3, 8, 9,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 9, 10,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 78, 0, 3, 10, 11,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 84, 0, 3, 11, 12,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 90, 0, 3, 12, 13,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 96, 0, 3, 13, 14,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 14, 15,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 15, 16,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 114, 0, 3, 16, 17,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 120, 0, 3, 17, 18,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 126, 0, 3, 18, 19,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 21, 24,
                                                                       60, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 24, 27,
                                                                       66, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 27, 30,
                                                                       72, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 162, 0, 3, 30, 33,
                                                                       78, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 33, 36,
                                                                       84, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 36, 39,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 39, 42,
                                                                       96, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 42, 45,
                                                                       102, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 45, 48,
                                                                       108, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 222, 0, 3, 48, 51,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 51, 54,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 242, 0, 3, 60, 66,
                                                                       132, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 257, 0, 3, 66, 72,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 272, 0, 3, 72, 78,
                                                                       152, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 287, 0, 3, 78, 84,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 302, 0, 3, 84, 90,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 317, 0, 3, 90, 96,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 332, 0, 3, 96,
                                                                       102, 192, 202, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 347, 0, 3, 102,
                                                                       108, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 362, 0, 3, 108,
                                                                       114, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 377, 0, 3, 114,
                                                                       120, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 392, 0, 3, 132,
                                                                       142, 242, 257, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 413, 0, 3, 142,
                                                                       152, 257, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 434, 0, 3, 152,
                                                                       162, 272, 287, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 455, 0, 3, 162,
                                                                       172, 287, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 476, 0, 3, 172,
                                                                       182, 302, 317, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 497, 0, 3, 182,
                                                                       192, 317, 332, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 192,
                                                                       202, 332, 347, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 539, 0, 3, 202,
                                                                       212, 347, 362, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 560, 0, 3, 212,
                                                                       222, 362, 377, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 581, 0, 3, 242,
                                                                       257, 392, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 609, 0, 3, 257,
                                                                       272, 413, 434, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 637, 0, 3, 272,
                                                                       287, 434, 455, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 665, 0, 3, 287,
                                                                       302, 455, 476, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 693, 0, 3, 302,
                                                                       317, 476, 497, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 721, 0, 3, 317,
                                                                       332, 497, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 749, 0, 3, 332,
                                                                       347, 518, 539, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 777, 0, 3, 347,
                                                                       362, 539, 560, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 805, 0, 3, 392,
                                                                       413, 581, 609, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 841, 0, 3, 413,
                                                                       434, 609, 637, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 877, 0, 3, 434,
                                                                       455, 637, 665, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 913, 0, 3, 455,
                                                                       476, 665, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 949, 0, 3, 476,
                                                                       497, 693, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 985, 0, 3, 497,
                                                                       518, 721, 749, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1021, 0, 3, 518,
                                                                       539, 749, 777, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1057, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1060, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1063, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1066, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1069, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1072, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1075, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1078, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1081, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1084, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1087, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1090, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1093, 3, 9, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1102, 3, 10, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1111, 3, 11, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1120, 3, 12, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1129, 3, 13, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1138, 3, 14, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1147, 3, 15, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1156, 3, 16, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1165, 3, 17, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1174, 3, 18, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1183, 3, 19, 57,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1192, 3, 27, 72,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1210, 3, 30, 78,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1228, 3, 33, 84,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1246, 3, 36, 90,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1264, 3, 39, 96,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1282, 3, 42, 102,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1300, 3, 45, 108,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1318, 3, 48, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1336, 3, 51, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1354, 3, 54, 126,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1372, 3, 72, 152,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1402, 3, 78, 162,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1432, 3, 84, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1462, 3, 90, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1492, 3, 96, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1522, 3, 102, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1552, 3, 108, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1582, 3, 114, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1612, 3, 120, 232,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1642, 3, 152, 272,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1687, 3, 162, 287,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1732, 3, 172, 302,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1777, 3, 182, 317,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1822, 3, 192, 332,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1867, 3, 202, 347,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1912, 3, 212, 362,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1957, 3, 222, 377,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2002, 3, 272, 434,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2065, 3, 287, 455,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2128, 3, 302, 476,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2191, 3, 317, 497,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2254, 3, 332, 518,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2317, 3, 347, 539,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2380, 3, 362, 560,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2443, 3, 434, 637,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2527, 3, 455, 665,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2611, 3, 476, 693,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2695, 3, 497, 721,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2779, 3, 518, 749,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2863, 3, 539, 777,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 2947, 3, 637, 877,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3055, 3, 665, 913,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3163, 3, 693, 949,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3271, 3, 721, 985,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3379, 3, 749,
                                                                       1021, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3487, 3, 7, 8,
                                                                       1057, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3493, 3, 8, 9,
                                                                       1060, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3499, 3, 9, 10,
                                                                       1063, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3505, 3, 10, 11,
                                                                       1066, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3511, 3, 11, 12,
                                                                       1069, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3517, 3, 12, 13,
                                                                       1072, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3523, 3, 13, 14,
                                                                       1075, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3529, 3, 14, 15,
                                                                       1078, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3535, 3, 15, 16,
                                                                       1081, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3541, 3, 16, 17,
                                                                       1084, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3547, 3, 17, 18,
                                                                       1087, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3553, 3, 18, 19,
                                                                       1090, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3559, 0, 3, 3487,
                                                                       1057, 3493, 1093, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3577, 0, 3, 3493,
                                                                       1060, 3499, 1102, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3595, 0, 3, 3499,
                                                                       1063, 3505, 1111, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3613, 0, 3, 3505,
                                                                       1066, 3511, 1120, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3631, 0, 3, 3511,
                                                                       1069, 3517, 1129, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3649, 0, 3, 3517,
                                                                       1072, 3523, 1138, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3667, 0, 3, 3523,
                                                                       1075, 3529, 1147, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3685, 0, 3, 3529,
                                                                       1078, 3535, 1156, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3703, 0, 3, 3535,
                                                                       1081, 3541, 1165, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3721, 0, 3, 3541,
                                                                       1084, 3547, 1174, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3739, 0, 3, 3547,
                                                                       1087, 3553, 1183, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3757, 0, 3, 3559,
                                                                       1093, 3577, 60, 66, 1192,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3793, 0, 3, 3577,
                                                                       1102, 3595, 66, 72, 1210,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3829, 0, 3, 3595,
                                                                       1111, 3613, 72, 78, 1228,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3865, 0, 3, 3613,
                                                                       1120, 3631, 78, 84, 1246,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3901, 0, 3, 3631,
                                                                       1129, 3649, 84, 90, 1264,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3937, 0, 3, 3649,
                                                                       1138, 3667, 90, 96, 1282,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3973, 0, 3, 3667,
                                                                       1147, 3685, 96, 102, 1300,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4009, 0, 3, 3685,
                                                                       1156, 3703, 102, 108,
                                                                       1318, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4045, 0, 3, 3703,
                                                                       1165, 3721, 108, 114,
                                                                       1336, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4081, 0, 3, 3721,
                                                                       1174, 3739, 114, 120,
                                                                       1354, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4117, 0, 3, 3757,
                                                                       1192, 3793, 132, 142,
                                                                       1372, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4177, 0, 3, 3793,
                                                                       1210, 3829, 142, 152,
                                                                       1402, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4237, 0, 3, 3829,
                                                                       1228, 3865, 152, 162,
                                                                       1432, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4297, 0, 3, 3865,
                                                                       1246, 3901, 162, 172,
                                                                       1462, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4357, 0, 3, 3901,
                                                                       1264, 3937, 172, 182,
                                                                       1492, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4417, 0, 3, 3937,
                                                                       1282, 3973, 182, 192,
                                                                       1522, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4477, 0, 3, 3973,
                                                                       1300, 4009, 192, 202,
                                                                       1552, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4537, 0, 3, 4009,
                                                                       1318, 4045, 202, 212,
                                                                       1582, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4597, 0, 3, 4045,
                                                                       1336, 4081, 212, 222,
                                                                       1612, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4657, 0, 3, 4117,
                                                                       1372, 4177, 242, 257,
                                                                       1642, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4747, 0, 3, 4177,
                                                                       1402, 4237, 257, 272,
                                                                       1687, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4837, 0, 3, 4237,
                                                                       1432, 4297, 272, 287,
                                                                       1732, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4927, 0, 3, 4297,
                                                                       1462, 4357, 287, 302,
                                                                       1777, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5017, 0, 3, 4357,
                                                                       1492, 4417, 302, 317,
                                                                       1822, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5107, 0, 3, 4417,
                                                                       1522, 4477, 317, 332,
                                                                       1867, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5197, 0, 3, 4477,
                                                                       1552, 4537, 332, 347,
                                                                       1912, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5287, 0, 3, 4537,
                                                                       1582, 4597, 347, 362,
                                                                       1957, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5377, 0, 3, 4657,
                                                                       1642, 4747, 392, 413,
                                                                       2002, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5503, 0, 3, 4747,
                                                                       1687, 4837, 413, 434,
                                                                       2065, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5629, 0, 3, 4837,
                                                                       1732, 4927, 434, 455,
                                                                       2128, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5755, 0, 3, 4927,
                                                                       1777, 5017, 455, 476,
                                                                       2191, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5881, 0, 3, 5017,
                                                                       1822, 5107, 476, 497,
                                                                       2254, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6007, 0, 3, 5107,
                                                                       1867, 5197, 497, 518,
                                                                       2317, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6133, 0, 3, 5197,
                                                                       1912, 5287, 518, 539,
                                                                       2380, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6259, 0, 3, 5377,
                                                                       2002, 5503, 581, 609,
                                                                       2443, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6427, 0, 3, 5503,
                                                                       2065, 5629, 609, 637,
                                                                       2527, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6595, 0, 3, 5629,
                                                                       2128, 5755, 637, 665,
                                                                       2611, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6763, 0, 3, 5755,
                                                                       2191, 5881, 665, 693,
                                                                       2695, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6931, 0, 3, 5881,
                                                                       2254, 6007, 693, 721,
                                                                       2779, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7099, 0, 3, 6007,
                                                                       2317, 6133, 721, 749,
                                                                       2863, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7267, 0, 3, 6259,
                                                                       2443, 6427, 805, 841,
                                                                       2947, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7483, 0, 3, 6427,
                                                                       2527, 6595, 841, 877,
                                                                       3055, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7699, 0, 3, 6595,
                                                                       2611, 6763, 877, 913,
                                                                       3163, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7915, 0, 3, 6763,
                                                                       2695, 6931, 913, 949,
                                                                       3271, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 8131, 0, 3, 6931,
                                                                       2779, 7099, 949, 985,
                                                                       3379, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8347, 3, 1057,
                                                                       1060, 3499, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8357, 3, 1060,
                                                                       1063, 3505, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8367, 3, 1063,
                                                                       1066, 3511, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8377, 3, 1066,
                                                                       1069, 3517, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8387, 3, 1069,
                                                                       1072, 3523, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8397, 3, 1072,
                                                                       1075, 3529, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8407, 3, 1075,
                                                                       1078, 3535, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8417, 3, 1078,
                                                                       1081, 3541, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8427, 3, 1081,
                                                                       1084, 3547, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8437, 3, 1084,
                                                                       1087, 3553, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8447, 0, 3, 8347,
                                                                       3499, 8357, 1093, 1102,
                                                                       3595, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8477, 0, 3, 8357,
                                                                       3505, 8367, 1102, 1111,
                                                                       3613, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8507, 0, 3, 8367,
                                                                       3511, 8377, 1111, 1120,
                                                                       3631, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8537, 0, 3, 8377,
                                                                       3517, 8387, 1120, 1129,
                                                                       3649, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8567, 0, 3, 8387,
                                                                       3523, 8397, 1129, 1138,
                                                                       3667, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8597, 0, 3, 8397,
                                                                       3529, 8407, 1138, 1147,
                                                                       3685, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8627, 0, 3, 8407,
                                                                       3535, 8417, 1147, 1156,
                                                                       3703, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8657, 0, 3, 8417,
                                                                       3541, 8427, 1156, 1165,
                                                                       3721, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8687, 0, 3, 8427,
                                                                       3547, 8437, 1165, 1174,
                                                                       3739, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8717, 0, 3, 8447,
                                                                       3595, 8477, 1192, 1210,
                                                                       3829, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8777, 0, 3, 8477,
                                                                       3613, 8507, 1210, 1228,
                                                                       3865, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8837, 0, 3, 8507,
                                                                       3631, 8537, 1228, 1246,
                                                                       3901, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8897, 0, 3, 8537,
                                                                       3649, 8567, 1246, 1264,
                                                                       3937, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8957, 0, 3, 8567,
                                                                       3667, 8597, 1264, 1282,
                                                                       3973, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9017, 0, 3, 8597,
                                                                       3685, 8627, 1282, 1300,
                                                                       4009, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9077, 0, 3, 8627,
                                                                       3703, 8657, 1300, 1318,
                                                                       4045, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9137, 0, 3, 8657,
                                                                       3721, 8687, 1318, 1336,
                                                                       4081, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9197, 0, 3, 8717,
                                                                       3829, 8777, 1372, 1402,
                                                                       4237, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9297, 0, 3, 8777,
                                                                       3865, 8837, 1402, 1432,
                                                                       4297, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9397, 0, 3, 8837,
                                                                       3901, 8897, 1432, 1462,
                                                                       4357, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9497, 0, 3, 8897,
                                                                       3937, 8957, 1462, 1492,
                                                                       4417, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9597, 0, 3, 8957,
                                                                       3973, 9017, 1492, 1522,
                                                                       4477, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9697, 0, 3, 9017,
                                                                       4009, 9077, 1522, 1552,
                                                                       4537, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9797, 0, 3, 9077,
                                                                       4045, 9137, 1552, 1582,
                                                                       4597, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9897, 0, 3, 9197,
                                                                       4237, 9297, 1642, 1687,
                                                                       4837, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10047, 0, 3, 9297,
                                                                       4297, 9397, 1687, 1732,
                                                                       4927, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10197, 0, 3, 9397,
                                                                       4357, 9497, 1732, 1777,
                                                                       5017, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10347, 0, 3, 9497,
                                                                       4417, 9597, 1777, 1822,
                                                                       5107, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10497, 0, 3, 9597,
                                                                       4477, 9697, 1822, 1867,
                                                                       5197, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10647, 0, 3, 9697,
                                                                       4537, 9797, 1867, 1912,
                                                                       5287, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 10797, 0, 3, 9897,
                                                                       4837, 10047, 2002, 2065,
                                                                       5629, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 11007, 0, 3,
                                                                       10047, 4927, 10197, 2065,
                                                                       2128, 5755, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 11217, 0, 3,
                                                                       10197, 5017, 10347, 2128,
                                                                       2191, 5881, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 11427, 0, 3,
                                                                       10347, 5107, 10497, 2191,
                                                                       2254, 6007, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 11637, 0, 3,
                                                                       10497, 5197, 10647, 2254,
                                                                       2317, 6133, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 11847, 0, 3,
                                                                       10797, 5629, 11007, 2443,
                                                                       2527, 6595, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 12127, 0, 3,
                                                                       11007, 5755, 11217, 2527,
                                                                       2611, 6763, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 12407, 0, 3,
                                                                       11217, 5881, 11427, 2611,
                                                                       2695, 6931, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 12687, 0, 3,
                                                                       11427, 6007, 11637, 2695,
                                                                       2779, 7099, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 12967, 0, 3,
                                                                       11847, 6595, 12127, 2947,
                                                                       3055, 7699, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 13327, 0, 3,
                                                                       12127, 6763, 12407, 3055,
                                                                       3163, 7915, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 13687, 0, 3,
                                                                       12407, 6931, 12687, 3163,
                                                                       3271, 8131, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14047, 3, 3487,
                                                                       3493, 8347, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14062, 3, 3493,
                                                                       3499, 8357, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14077, 3, 3499,
                                                                       3505, 8367, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14092, 3, 3505,
                                                                       3511, 8377, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14107, 3, 3511,
                                                                       3517, 8387, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14122, 3, 3517,
                                                                       3523, 8397, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14137, 3, 3523,
                                                                       3529, 8407, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14152, 3, 3529,
                                                                       3535, 8417, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14167, 3, 3535,
                                                                       3541, 8427, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14182, 3, 3541,
                                                                       3547, 8437, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14197, 0, 3,
                                                                       14047, 8347, 14062, 3559,
                                                                       3577, 8447, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14242, 0, 3,
                                                                       14062, 8357, 14077, 3577,
                                                                       3595, 8477, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14287, 0, 3,
                                                                       14077, 8367, 14092, 3595,
                                                                       3613, 8507, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14332, 0, 3,
                                                                       14092, 8377, 14107, 3613,
                                                                       3631, 8537, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14377, 0, 3,
                                                                       14107, 8387, 14122, 3631,
                                                                       3649, 8567, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14422, 0, 3,
                                                                       14122, 8397, 14137, 3649,
                                                                       3667, 8597, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14467, 0, 3,
                                                                       14137, 8407, 14152, 3667,
                                                                       3685, 8627, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14512, 0, 3,
                                                                       14152, 8417, 14167, 3685,
                                                                       3703, 8657, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14557, 0, 3,
                                                                       14167, 8427, 14182, 3703,
                                                                       3721, 8687, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14602, 0, 3,
                                                                       14197, 8447, 14242, 3757,
                                                                       3793, 8717, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14692, 0, 3,
                                                                       14242, 8477, 14287, 3793,
                                                                       3829, 8777, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14782, 0, 3,
                                                                       14287, 8507, 14332, 3829,
                                                                       3865, 8837, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14872, 0, 3,
                                                                       14332, 8537, 14377, 3865,
                                                                       3901, 8897, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14962, 0, 3,
                                                                       14377, 8567, 14422, 3901,
                                                                       3937, 8957, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 15052, 0, 3,
                                                                       14422, 8597, 14467, 3937,
                                                                       3973, 9017, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 15142, 0, 3,
                                                                       14467, 8627, 14512, 3973,
                                                                       4009, 9077, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 15232, 0, 3,
                                                                       14512, 8657, 14557, 4009,
                                                                       4045, 9137, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15322, 0, 3,
                                                                       14602, 8717, 14692, 4117,
                                                                       4177, 9197, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15472, 0, 3,
                                                                       14692, 8777, 14782, 4177,
                                                                       4237, 9297, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15622, 0, 3,
                                                                       14782, 8837, 14872, 4237,
                                                                       4297, 9397, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15772, 0, 3,
                                                                       14872, 8897, 14962, 4297,
                                                                       4357, 9497, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15922, 0, 3,
                                                                       14962, 8957, 15052, 4357,
                                                                       4417, 9597, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 16072, 0, 3,
                                                                       15052, 9017, 15142, 4417,
                                                                       4477, 9697, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 16222, 0, 3,
                                                                       15142, 9077, 15232, 4477,
                                                                       4537, 9797, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 16372, 0, 3,
                                                                       15322, 9197, 15472, 4657,
                                                                       4747, 9897, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 16597, 0, 3,
                                                                       15472, 9297, 15622, 4747,
                                                                       4837, 10047, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 16822, 0, 3,
                                                                       15622, 9397, 15772, 4837,
                                                                       4927, 10197, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 17047, 0, 3,
                                                                       15772, 9497, 15922, 4927,
                                                                       5017, 10347, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 17272, 0, 3,
                                                                       15922, 9597, 16072, 5017,
                                                                       5107, 10497, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 17497, 0, 3,
                                                                       16072, 9697, 16222, 5107,
                                                                       5197, 10647, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 17722, 0, 3,
                                                                       16372, 9897, 16597, 5377,
                                                                       5503, 10797, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 18037, 0, 3,
                                                                       16597, 10047, 16822, 5503,
                                                                       5629, 11007, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 18352, 0, 3,
                                                                       16822, 10197, 17047, 5629,
                                                                       5755, 11217, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 18667, 0, 3,
                                                                       17047, 10347, 17272, 5755,
                                                                       5881, 11427, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 18982, 0, 3,
                                                                       17272, 10497, 17497, 5881,
                                                                       6007, 11637, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 19297, 0, 3,
                                                                       17722, 10797, 18037, 6259,
                                                                       6427, 11847, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 19717, 0, 3,
                                                                       18037, 11007, 18352, 6427,
                                                                       6595, 12127, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 20137, 0, 3,
                                                                       18352, 11217, 18667, 6595,
                                                                       6763, 12407, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 20557, 0, 3,
                                                                       18667, 11427, 18982, 6763,
                                                                       6931, 12687, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 20977, 0, 3,
                                                                       19297, 11847, 19717, 7267,
                                                                       7483, 12967, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 21517, 0, 3,
                                                                       19717, 12127, 20137, 7483,
                                                                       7699, 13327, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 22057, 0, 3,
                                                                       20137, 12407, 20557, 7699,
                                                                       7915, 13687, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22597, 3, 8347,
                                                                       8357, 14077, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22618, 3, 8357,
                                                                       8367, 14092, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22639, 3, 8367,
                                                                       8377, 14107, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22660, 3, 8377,
                                                                       8387, 14122, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22681, 3, 8387,
                                                                       8397, 14137, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22702, 3, 8397,
                                                                       8407, 14152, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22723, 3, 8407,
                                                                       8417, 14167, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22744, 3, 8417,
                                                                       8427, 14182, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 22765, 0, 3,
                                                                       22597, 14077, 22618, 8447,
                                                                       8477, 14287, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 22828, 0, 3,
                                                                       22618, 14092, 22639, 8477,
                                                                       8507, 14332, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 22891, 0, 3,
                                                                       22639, 14107, 22660, 8507,
                                                                       8537, 14377, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 22954, 0, 3,
                                                                       22660, 14122, 22681, 8537,
                                                                       8567, 14422, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23017, 0, 3,
                                                                       22681, 14137, 22702, 8567,
                                                                       8597, 14467, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23080, 0, 3,
                                                                       22702, 14152, 22723, 8597,
                                                                       8627, 14512, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23143, 0, 3,
                                                                       22723, 14167, 22744, 8627,
                                                                       8657, 14557, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23206, 0, 3,
                                                                       22765, 14287, 22828, 8717,
                                                                       8777, 14782, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23332, 0, 3,
                                                                       22828, 14332, 22891, 8777,
                                                                       8837, 14872, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23458, 0, 3,
                                                                       22891, 14377, 22954, 8837,
                                                                       8897, 14962, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23584, 0, 3,
                                                                       22954, 14422, 23017, 8897,
                                                                       8957, 15052, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23710, 0, 3,
                                                                       23017, 14467, 23080, 8957,
                                                                       9017, 15142, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23836, 0, 3,
                                                                       23080, 14512, 23143, 9017,
                                                                       9077, 15232, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 23962, 0, 3,
                                                                       23206, 14782, 23332, 9197,
                                                                       9297, 15622, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 24172, 0, 3,
                                                                       23332, 14872, 23458, 9297,
                                                                       9397, 15772, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 24382, 0, 3,
                                                                       23458, 14962, 23584, 9397,
                                                                       9497, 15922, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 24592, 0, 3,
                                                                       23584, 15052, 23710, 9497,
                                                                       9597, 16072, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 24802, 0, 3,
                                                                       23710, 15142, 23836, 9597,
                                                                       9697, 16222, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 25012, 0, 3,
                                                                       23962, 15622, 24172, 9897,
                                                                       10047, 16822, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 25327, 0, 3,
                                                                       24172, 15772, 24382,
                                                                       10047, 10197, 17047,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 25642, 0, 3,
                                                                       24382, 15922, 24592,
                                                                       10197, 10347, 17272,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 25957, 0, 3,
                                                                       24592, 16072, 24802,
                                                                       10347, 10497, 17497,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 26272, 0, 3,
                                                                       25012, 16822, 25327,
                                                                       10797, 11007, 18352,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 26713, 0, 3,
                                                                       25327, 17047, 25642,
                                                                       11007, 11217, 18667,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 27154, 0, 3,
                                                                       25642, 17272, 25957,
                                                                       11217, 11427, 18982,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 27595, 0, 3,
                                                                       26272, 18352, 26713,
                                                                       11847, 12127, 20137,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 28183, 0, 3,
                                                                       26713, 18667, 27154,
                                                                       12127, 12407, 20557,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 28771, 0, 3,
                                                                       27595, 20137, 28183,
                                                                       12967, 13327, 22057,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29527, 3, 14047,
                                                                       14062, 22597, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29555, 3, 14062,
                                                                       14077, 22618, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29583, 3, 14077,
                                                                       14092, 22639, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29611, 3, 14092,
                                                                       14107, 22660, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29639, 3, 14107,
                                                                       14122, 22681, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29667, 3, 14122,
                                                                       14137, 22702, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29695, 3, 14137,
                                                                       14152, 22723, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29723, 3, 14152,
                                                                       14167, 22744, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29751, 0, 3,
                                                                       29527, 22597, 29555,
                                                                       14197, 14242, 22765,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29835, 0, 3,
                                                                       29555, 22618, 29583,
                                                                       14242, 14287, 22828,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29919, 0, 3,
                                                                       29583, 22639, 29611,
                                                                       14287, 14332, 22891,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 30003, 0, 3,
                                                                       29611, 22660, 29639,
                                                                       14332, 14377, 22954,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 30087, 0, 3,
                                                                       29639, 22681, 29667,
                                                                       14377, 14422, 23017,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 30171, 0, 3,
                                                                       29667, 22702, 29695,
                                                                       14422, 14467, 23080,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 30255, 0, 3,
                                                                       29695, 22723, 29723,
                                                                       14467, 14512, 23143,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 30339, 0, 3,
                                                                       29751, 22765, 29835,
                                                                       14602, 14692, 23206,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 30507, 0, 3,
                                                                       29835, 22828, 29919,
                                                                       14692, 14782, 23332,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 30675, 0, 3,
                                                                       29919, 22891, 30003,
                                                                       14782, 14872, 23458,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 30843, 0, 3,
                                                                       30003, 22954, 30087,
                                                                       14872, 14962, 23584,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 31011, 0, 3,
                                                                       30087, 23017, 30171,
                                                                       14962, 15052, 23710,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 31179, 0, 3,
                                                                       30171, 23080, 30255,
                                                                       15052, 15142, 23836,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 31347, 0, 3,
                                                                       30339, 23206, 30507,
                                                                       15322, 15472, 23962,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 31627, 0, 3,
                                                                       30507, 23332, 30675,
                                                                       15472, 15622, 24172,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 31907, 0, 3,
                                                                       30675, 23458, 30843,
                                                                       15622, 15772, 24382,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 32187, 0, 3,
                                                                       30843, 23584, 31011,
                                                                       15772, 15922, 24592,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 32467, 0, 3,
                                                                       31011, 23710, 31179,
                                                                       15922, 16072, 24802,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 32747, 0, 3,
                                                                       31347, 23962, 31627,
                                                                       16372, 16597, 25012,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 33167, 0, 3,
                                                                       31627, 24172, 31907,
                                                                       16597, 16822, 25327,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 33587, 0, 3,
                                                                       31907, 24382, 32187,
                                                                       16822, 17047, 25642,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 34007, 0, 3,
                                                                       32187, 24592, 32467,
                                                                       17047, 17272, 25957,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 34427, 0, 3,
                                                                       32747, 25012, 33167,
                                                                       17722, 18037, 26272,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 35015, 0, 3,
                                                                       33167, 25327, 33587,
                                                                       18037, 18352, 26713,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 35603, 0, 3,
                                                                       33587, 25642, 34007,
                                                                       18352, 18667, 27154,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 36191, 0, 3,
                                                                       34427, 26272, 35015,
                                                                       19297, 19717, 27595,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 36975, 0, 3,
                                                                       35015, 26713, 35603,
                                                                       19717, 20137, 28183,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 37759, 0, 3,
                                                                       36191, 27595, 36975,
                                                                       20977, 21517, 28771,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 38767, 34427, 588, ncols);

                    simdfunc::contract_primitives(buffer, 39628, 36191, 784, ncols);

                    simdfunc::contract_primitives(buffer, 40776, 37759, 1008, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 39355, 38767, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 40412, 39628, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 41784, 40776, 36, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 42252, 39355, 40412, 13,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 43071, 40412, 41784, 13,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 44163, 42252, 43071, 13,
                                             nmax);

        simdtrf::transform_d_inner(buffer, 45801, 44163, 21, 13, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 45801, 65, nmax);
    }

    for (size_t m = 0; m < 715; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
