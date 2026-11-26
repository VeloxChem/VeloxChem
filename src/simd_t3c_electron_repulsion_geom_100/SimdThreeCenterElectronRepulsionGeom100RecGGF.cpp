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


#include "SimdThreeCenterElectronRepulsionGeom100RecGGF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferGD.hpp"
#include "SimdTransferGF.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferGeom100XGD.hpp"
#include "SimdTransferGeom100XGF.hpp"
#include "SimdTransferGeom100XGG.hpp"
#include "SimdTransferGeom100XGP.hpp"
#include "SimdTransferGeom100XHD.hpp"
#include "SimdTransferGeom100XHF.hpp"
#include "SimdTransferGeom100XHP.hpp"
#include "SimdTransferGeom100XID.hpp"
#include "SimdTransferGeom100XIP.hpp"
#include "SimdTransferGeom100XKP.hpp"
#include "SimdTransferGeom100YGD.hpp"
#include "SimdTransferGeom100YGF.hpp"
#include "SimdTransferGeom100YGG.hpp"
#include "SimdTransferGeom100YGP.hpp"
#include "SimdTransferGeom100YHD.hpp"
#include "SimdTransferGeom100YHF.hpp"
#include "SimdTransferGeom100YHP.hpp"
#include "SimdTransferGeom100YID.hpp"
#include "SimdTransferGeom100YIP.hpp"
#include "SimdTransferGeom100YKP.hpp"
#include "SimdTransferGeom100ZGD.hpp"
#include "SimdTransferGeom100ZGF.hpp"
#include "SimdTransferGeom100ZGG.hpp"
#include "SimdTransferGeom100ZGP.hpp"
#include "SimdTransferGeom100ZHD.hpp"
#include "SimdTransferGeom100ZHF.hpp"
#include "SimdTransferGeom100ZHP.hpp"
#include "SimdTransferGeom100ZID.hpp"
#include "SimdTransferGeom100ZIP.hpp"
#include "SimdTransferGeom100ZKP.hpp"
#include "SimdTransferHD.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_ggf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_ggf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 61012, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 61012, 20417, 8780, dimensions);

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 6, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 19, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 22, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 7, 8,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 8, 9,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 9, 10,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 10, 11,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 76, 0, 3, 11, 12,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 12, 13,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 88, 0, 3, 13, 14,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 94, 0, 3, 14, 15,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 100, 0, 3, 15, 16,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 106, 0, 3, 16, 17,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 19, 22,
                                                                       52, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 22, 25,
                                                                       58, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 25, 28,
                                                                       64, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 28, 31,
                                                                       70, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 31, 34,
                                                                       76, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 162, 0, 3, 34, 37,
                                                                       82, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 37, 40,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 40, 43,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 43, 46,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 52, 58,
                                                                       112, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 217, 0, 3, 58, 64,
                                                                       122, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 64, 70,
                                                                       132, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 247, 0, 3, 70, 76,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 262, 0, 3, 76, 82,
                                                                       152, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 277, 0, 3, 82, 88,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 292, 0, 3, 88, 94,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 307, 0, 3, 94,
                                                                       100, 182, 192, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 322, 0, 3, 112,
                                                                       122, 202, 217, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 343, 0, 3, 122,
                                                                       132, 217, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 364, 0, 3, 132,
                                                                       142, 232, 247, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 385, 0, 3, 142,
                                                                       152, 247, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 406, 0, 3, 152,
                                                                       162, 262, 277, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 427, 0, 3, 162,
                                                                       172, 277, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 172,
                                                                       182, 292, 307, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 469, 0, 3, 202,
                                                                       217, 322, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 497, 0, 3, 217,
                                                                       232, 343, 364, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 525, 0, 3, 232,
                                                                       247, 364, 385, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 553, 0, 3, 247,
                                                                       262, 385, 406, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 581, 0, 3, 262,
                                                                       277, 406, 427, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 609, 0, 3, 277,
                                                                       292, 427, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 637, 0, 3, 322,
                                                                       343, 469, 497, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 673, 0, 3, 343,
                                                                       364, 497, 525, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 709, 0, 3, 364,
                                                                       385, 525, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 745, 0, 3, 385,
                                                                       406, 553, 581, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 781, 0, 3, 406,
                                                                       427, 581, 609, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 817, 0, 3, 469,
                                                                       497, 637, 673, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 862, 0, 3, 497,
                                                                       525, 673, 709, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 907, 0, 3, 525,
                                                                       553, 709, 745, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 952, 0, 3, 553,
                                                                       581, 745, 781, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 997, 0, 3, 637,
                                                                       673, 817, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1052, 0, 3, 673,
                                                                       709, 862, 907, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1107, 0, 3, 709,
                                                                       745, 907, 952, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1162, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1165, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1168, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1171, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1174, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1177, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1180, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1183, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1186, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1189, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1192, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1195, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1198, 3, 7, 19,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1207, 3, 8, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1216, 3, 9, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1225, 3, 10, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1234, 3, 11, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1243, 3, 12, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1252, 3, 13, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1261, 3, 14, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1270, 3, 15, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1279, 3, 16, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1288, 3, 17, 49,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1297, 3, 19, 52,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1315, 3, 22, 58,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1333, 3, 25, 64,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1351, 3, 28, 70,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1369, 3, 31, 76,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1387, 3, 34, 82,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1405, 3, 37, 88,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1423, 3, 40, 94,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1441, 3, 43, 100,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1459, 3, 46, 106,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1477, 3, 52, 112,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1507, 3, 58, 122,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1537, 3, 64, 132,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1567, 3, 70, 142,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1597, 3, 76, 152,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1627, 3, 82, 162,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1657, 3, 88, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1687, 3, 94, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1717, 3, 100, 192,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1747, 3, 112, 202,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1792, 3, 122, 217,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1837, 3, 132, 232,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1882, 3, 142, 247,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1927, 3, 152, 262,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1972, 3, 162, 277,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2017, 3, 172, 292,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2062, 3, 182, 307,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2107, 3, 202, 322,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2170, 3, 217, 343,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2233, 3, 232, 364,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2296, 3, 247, 385,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2359, 3, 262, 406,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2422, 3, 277, 427,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2485, 3, 292, 448,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2548, 3, 322, 469,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2632, 3, 343, 497,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2716, 3, 364, 525,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2800, 3, 385, 553,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2884, 3, 406, 581,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2968, 3, 427, 609,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3052, 3, 469, 637,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3160, 3, 497, 673,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3268, 3, 525, 709,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3376, 3, 553, 745,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3484, 3, 581, 781,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 3592, 3, 637, 817,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 3727, 3, 673, 862,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 3862, 3, 709, 907,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 3997, 3, 745, 952,
                                                                       ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 4132, 3, 817, 997,
                                                                       ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 4297, 3, 862,
                                                                       1052, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 4462, 3, 907,
                                                                       1107, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4627, 3, 7, 8,
                                                                       1168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4633, 3, 8, 9,
                                                                       1171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4639, 3, 9, 10,
                                                                       1174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4645, 3, 10, 11,
                                                                       1177, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4651, 3, 11, 12,
                                                                       1180, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4657, 3, 12, 13,
                                                                       1183, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4663, 3, 13, 14,
                                                                       1186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4669, 3, 14, 15,
                                                                       1189, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4675, 3, 15, 16,
                                                                       1192, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4681, 3, 16, 17,
                                                                       1195, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4687, 0, 3, 4627,
                                                                       1168, 4633, 1216, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4705, 0, 3, 4633,
                                                                       1171, 4639, 1225, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4723, 0, 3, 4639,
                                                                       1174, 4645, 1234, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4741, 0, 3, 4645,
                                                                       1177, 4651, 1243, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4759, 0, 3, 4651,
                                                                       1180, 4657, 1252, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4777, 0, 3, 4657,
                                                                       1183, 4663, 1261, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4795, 0, 3, 4663,
                                                                       1186, 4669, 1270, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4813, 0, 3, 4669,
                                                                       1189, 4675, 1279, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4831, 0, 3, 4675,
                                                                       1192, 4681, 1288, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4849, 0, 3, 4687,
                                                                       1216, 4705, 52, 58, 1333,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4885, 0, 3, 4705,
                                                                       1225, 4723, 58, 64, 1351,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4921, 0, 3, 4723,
                                                                       1234, 4741, 64, 70, 1369,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4957, 0, 3, 4741,
                                                                       1243, 4759, 70, 76, 1387,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4993, 0, 3, 4759,
                                                                       1252, 4777, 76, 82, 1405,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5029, 0, 3, 4777,
                                                                       1261, 4795, 82, 88, 1423,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5065, 0, 3, 4795,
                                                                       1270, 4813, 88, 94, 1441,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5101, 0, 3, 4813,
                                                                       1279, 4831, 94, 100, 1459,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5137, 0, 3, 4849,
                                                                       1333, 4885, 112, 122,
                                                                       1537, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5197, 0, 3, 4885,
                                                                       1351, 4921, 122, 132,
                                                                       1567, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5257, 0, 3, 4921,
                                                                       1369, 4957, 132, 142,
                                                                       1597, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5317, 0, 3, 4957,
                                                                       1387, 4993, 142, 152,
                                                                       1627, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5377, 0, 3, 4993,
                                                                       1405, 5029, 152, 162,
                                                                       1657, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5437, 0, 3, 5029,
                                                                       1423, 5065, 162, 172,
                                                                       1687, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5497, 0, 3, 5065,
                                                                       1441, 5101, 172, 182,
                                                                       1717, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5557, 0, 3, 5137,
                                                                       1537, 5197, 202, 217,
                                                                       1837, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5647, 0, 3, 5197,
                                                                       1567, 5257, 217, 232,
                                                                       1882, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5737, 0, 3, 5257,
                                                                       1597, 5317, 232, 247,
                                                                       1927, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5827, 0, 3, 5317,
                                                                       1627, 5377, 247, 262,
                                                                       1972, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5917, 0, 3, 5377,
                                                                       1657, 5437, 262, 277,
                                                                       2017, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6007, 0, 3, 5437,
                                                                       1687, 5497, 277, 292,
                                                                       2062, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6097, 0, 3, 5557,
                                                                       1837, 5647, 322, 343,
                                                                       2233, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6223, 0, 3, 5647,
                                                                       1882, 5737, 343, 364,
                                                                       2296, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6349, 0, 3, 5737,
                                                                       1927, 5827, 364, 385,
                                                                       2359, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6475, 0, 3, 5827,
                                                                       1972, 5917, 385, 406,
                                                                       2422, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6601, 0, 3, 5917,
                                                                       2017, 6007, 406, 427,
                                                                       2485, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6727, 0, 3, 6097,
                                                                       2233, 6223, 469, 497,
                                                                       2716, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6895, 0, 3, 6223,
                                                                       2296, 6349, 497, 525,
                                                                       2800, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7063, 0, 3, 6349,
                                                                       2359, 6475, 525, 553,
                                                                       2884, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7231, 0, 3, 6475,
                                                                       2422, 6601, 553, 581,
                                                                       2968, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7399, 0, 3, 6727,
                                                                       2716, 6895, 637, 673,
                                                                       3268, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7615, 0, 3, 6895,
                                                                       2800, 7063, 673, 709,
                                                                       3376, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 7831, 0, 3, 7063,
                                                                       2884, 7231, 709, 745,
                                                                       3484, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 8047, 0, 3, 7399,
                                                                       3268, 7615, 817, 862,
                                                                       3862, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 8317, 0, 3, 7615,
                                                                       3376, 7831, 862, 907,
                                                                       3997, ncols, gamma, p,
                                                                       q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 8587, 0, 3, 8047,
                                                                       3862, 8317, 997, 1052,
                                                                       4462, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8917, 3, 1162,
                                                                       1165, 4627, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8927, 3, 1165,
                                                                       1168, 4633, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8937, 3, 1168,
                                                                       1171, 4639, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8947, 3, 1171,
                                                                       1174, 4645, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8957, 3, 1174,
                                                                       1177, 4651, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8967, 3, 1177,
                                                                       1180, 4657, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8977, 3, 1180,
                                                                       1183, 4663, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8987, 3, 1183,
                                                                       1186, 4669, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8997, 3, 1186,
                                                                       1189, 4675, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9007, 3, 1189,
                                                                       1192, 4681, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9017, 0, 3, 8917,
                                                                       4627, 8927, 1198, 1207,
                                                                       4687, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9047, 0, 3, 8927,
                                                                       4633, 8937, 1207, 1216,
                                                                       4705, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9077, 0, 3, 8937,
                                                                       4639, 8947, 1216, 1225,
                                                                       4723, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9107, 0, 3, 8947,
                                                                       4645, 8957, 1225, 1234,
                                                                       4741, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9137, 0, 3, 8957,
                                                                       4651, 8967, 1234, 1243,
                                                                       4759, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9167, 0, 3, 8967,
                                                                       4657, 8977, 1243, 1252,
                                                                       4777, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9197, 0, 3, 8977,
                                                                       4663, 8987, 1252, 1261,
                                                                       4795, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9227, 0, 3, 8987,
                                                                       4669, 8997, 1261, 1270,
                                                                       4813, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9257, 0, 3, 8997,
                                                                       4675, 9007, 1270, 1279,
                                                                       4831, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9287, 0, 3, 9017,
                                                                       4687, 9047, 1297, 1315,
                                                                       4849, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9347, 0, 3, 9047,
                                                                       4705, 9077, 1315, 1333,
                                                                       4885, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9407, 0, 3, 9077,
                                                                       4723, 9107, 1333, 1351,
                                                                       4921, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9467, 0, 3, 9107,
                                                                       4741, 9137, 1351, 1369,
                                                                       4957, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9527, 0, 3, 9137,
                                                                       4759, 9167, 1369, 1387,
                                                                       4993, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9587, 0, 3, 9167,
                                                                       4777, 9197, 1387, 1405,
                                                                       5029, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9647, 0, 3, 9197,
                                                                       4795, 9227, 1405, 1423,
                                                                       5065, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9707, 0, 3, 9227,
                                                                       4813, 9257, 1423, 1441,
                                                                       5101, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9767, 0, 3, 9287,
                                                                       4849, 9347, 1477, 1507,
                                                                       5137, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9867, 0, 3, 9347,
                                                                       4885, 9407, 1507, 1537,
                                                                       5197, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9967, 0, 3, 9407,
                                                                       4921, 9467, 1537, 1567,
                                                                       5257, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10067, 0, 3, 9467,
                                                                       4957, 9527, 1567, 1597,
                                                                       5317, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10167, 0, 3, 9527,
                                                                       4993, 9587, 1597, 1627,
                                                                       5377, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10267, 0, 3, 9587,
                                                                       5029, 9647, 1627, 1657,
                                                                       5437, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10367, 0, 3, 9647,
                                                                       5065, 9707, 1657, 1687,
                                                                       5497, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10467, 0, 3, 9767,
                                                                       5137, 9867, 1747, 1792,
                                                                       5557, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10617, 0, 3, 9867,
                                                                       5197, 9967, 1792, 1837,
                                                                       5647, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10767, 0, 3, 9967,
                                                                       5257, 10067, 1837, 1882,
                                                                       5737, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 10917, 0, 3,
                                                                       10067, 5317, 10167, 1882,
                                                                       1927, 5827, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11067, 0, 3,
                                                                       10167, 5377, 10267, 1927,
                                                                       1972, 5917, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11217, 0, 3,
                                                                       10267, 5437, 10367, 1972,
                                                                       2017, 6007, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 11367, 0, 3,
                                                                       10467, 5557, 10617, 2107,
                                                                       2170, 6097, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 11577, 0, 3,
                                                                       10617, 5647, 10767, 2170,
                                                                       2233, 6223, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 11787, 0, 3,
                                                                       10767, 5737, 10917, 2233,
                                                                       2296, 6349, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 11997, 0, 3,
                                                                       10917, 5827, 11067, 2296,
                                                                       2359, 6475, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 12207, 0, 3,
                                                                       11067, 5917, 11217, 2359,
                                                                       2422, 6601, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 12417, 0, 3,
                                                                       11367, 6097, 11577, 2548,
                                                                       2632, 6727, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 12697, 0, 3,
                                                                       11577, 6223, 11787, 2632,
                                                                       2716, 6895, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 12977, 0, 3,
                                                                       11787, 6349, 11997, 2716,
                                                                       2800, 7063, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 13257, 0, 3,
                                                                       11997, 6475, 12207, 2800,
                                                                       2884, 7231, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 13537, 0, 3,
                                                                       12417, 6727, 12697, 3052,
                                                                       3160, 7399, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 13897, 0, 3,
                                                                       12697, 6895, 12977, 3160,
                                                                       3268, 7615, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 14257, 0, 3,
                                                                       12977, 7063, 13257, 3268,
                                                                       3376, 7831, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 14617, 0, 3,
                                                                       13537, 7399, 13897, 3592,
                                                                       3727, 8047, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 15067, 0, 3,
                                                                       13897, 7615, 14257, 3727,
                                                                       3862, 8317, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 15517, 0, 3,
                                                                       14617, 8047, 15067, 4132,
                                                                       4297, 8587, ncols, gamma,
                                                                       p, q);

                    simdgeo::geom_g_x(buffer, 16067, 9767, 11367, 1, 10, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 16217, 9767, 11367, 1, 10, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 16367, 9767, 11367, 1, 10, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 16517, 10467, 12417, 1, 10, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 16727, 10467, 12417, 1, 10, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 16937, 10467, 12417, 1, 10, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 17147, 11367, 13537, 1, 10, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 17427, 11367, 13537, 1, 10, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 17707, 11367, 13537, 1, 10, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 17987, 12417, 14617, 1, 10, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 18347, 12417, 14617, 1, 10, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 18707, 12417, 14617, 1, 10, ncols, alpha);

                    simdgeo::geom_l_x(buffer, 19067, 13537, 15517, 1, 10, ncols, alpha);

                    simdgeo::geom_l_y(buffer, 19517, 13537, 15517, 1, 10, ncols, alpha);

                    simdgeo::geom_l_z(buffer, 19967, 13537, 15517, 1, 10, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 20417, 16067, 150, ncols);

                    simdfunc::contract_primitives(buffer, 20672, 16217, 150, ncols);

                    simdfunc::contract_primitives(buffer, 20927, 16367, 150, ncols);

                    simdfunc::contract_primitives(buffer, 21182, 10467, 150, ncols);

                    simdfunc::contract_primitives(buffer, 21437, 16517, 210, ncols);

                    simdfunc::contract_primitives(buffer, 21794, 16727, 210, ncols);

                    simdfunc::contract_primitives(buffer, 22151, 16937, 210, ncols);

                    simdfunc::contract_primitives(buffer, 22508, 11367, 210, ncols);

                    simdfunc::contract_primitives(buffer, 22865, 17147, 280, ncols);

                    simdfunc::contract_primitives(buffer, 23341, 17427, 280, ncols);

                    simdfunc::contract_primitives(buffer, 23817, 17707, 280, ncols);

                    simdfunc::contract_primitives(buffer, 24293, 12417, 280, ncols);

                    simdfunc::contract_primitives(buffer, 24769, 17987, 360, ncols);

                    simdfunc::contract_primitives(buffer, 25381, 18347, 360, ncols);

                    simdfunc::contract_primitives(buffer, 25993, 18707, 360, ncols);

                    simdfunc::contract_primitives(buffer, 26605, 13537, 360, ncols);

                    simdfunc::contract_primitives(buffer, 27217, 19067, 450, ncols);

                    simdfunc::contract_primitives(buffer, 27982, 19517, 450, ncols);

                    simdfunc::contract_primitives(buffer, 28747, 19967, 450, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 20567, 20417, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 20822, 20672, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 21077, 20927, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 21332, 21182, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 21647, 21437, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 22004, 21794, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 22361, 22151, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 22718, 22508, 21, 1, nmax);

        simdtrf::transform_f_inner(buffer, 23145, 22865, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 23621, 23341, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 24097, 23817, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 24573, 24293, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 25129, 24769, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 25741, 25381, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 26353, 25993, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 26965, 26605, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 27667, 27217, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 28432, 27982, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 29197, 28747, 45, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 29512, 20567, 21332,
                                                       21647, 7, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 29827, 20822, 21332,
                                                       22004, 7, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 30142, 21077, 21332,
                                                       22361, 7, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 30457, 21332, 22718, 7, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 30772, 21647, 22718,
                                                       23145, 7, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 31213, 22004, 22718,
                                                       23621, 7, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 31654, 22361, 22718,
                                                       24097, 7, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 32095, 22718, 24573, 7, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 32536, 23145, 24573,
                                                       25129, 7, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 33124, 23621, 24573,
                                                       25741, 7, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 33712, 24097, 24573,
                                                       26353, 7, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 34300, 24573, 26965, 7, nmax);

        simdtrf::compute_hrr_geom_100x_kp_out_of_first(buffer, coordinates, 34888, 25129, 26965,
                                                       27667, 7, nmax);

        simdtrf::compute_hrr_geom_100y_kp_out_of_first(buffer, coordinates, 35644, 25741, 26965,
                                                       28432, 7, nmax);

        simdtrf::compute_hrr_geom_100z_kp_out_of_first(buffer, coordinates, 36400, 26353, 26965,
                                                       29197, 7, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 37156, 29512, 30457,
                                                       30772, 7, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 37786, 29827, 30457,
                                                       31213, 7, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 38416, 30142, 30457,
                                                       31654, 7, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 39046, 30457, 32095, 7, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 39676, 30772, 32095,
                                                       32536, 7, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 40558, 31213, 32095,
                                                       33124, 7, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 41440, 31654, 32095,
                                                       33712, 7, nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 42322, 32095, 34300, 7, nmax);

        simdtrf::compute_hrr_geom_100x_id_out_of_first(buffer, coordinates, 43204, 32536, 34300,
                                                       34888, 7, nmax);

        simdtrf::compute_hrr_geom_100y_id_out_of_first(buffer, coordinates, 44380, 33124, 34300,
                                                       35644, 7, nmax);

        simdtrf::compute_hrr_geom_100z_id_out_of_first(buffer, coordinates, 45556, 33712, 34300,
                                                       36400, 7, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 46732, 37156, 39046,
                                                       39676, 7, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 47782, 37786, 39046,
                                                       40558, 7, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 48832, 38416, 39046,
                                                       41440, 7, nmax);

        simdtrf::compute_hrr_gf_out_of_first(buffer, coordinates, 49882, 39046, 42322, 7, nmax);

        simdtrf::compute_hrr_geom_100x_hf_out_of_first(buffer, coordinates, 50932, 39676, 42322,
                                                       43204, 7, nmax);

        simdtrf::compute_hrr_geom_100y_hf_out_of_first(buffer, coordinates, 52402, 40558, 42322,
                                                       44380, 7, nmax);

        simdtrf::compute_hrr_geom_100z_hf_out_of_first(buffer, coordinates, 53872, 41440, 42322,
                                                       45556, 7, nmax);

        simdtrf::compute_hrr_geom_100x_gg_out_of_first(buffer, coordinates, 55342, 46732, 49882,
                                                       50932, 7, nmax);

        simdtrf::compute_hrr_geom_100y_gg_out_of_first(buffer, coordinates, 56917, 47782, 49882,
                                                       52402, 7, nmax);

        simdtrf::compute_hrr_geom_100z_gg_out_of_first(buffer, coordinates, 58492, 48832, 49882,
                                                       53872, 7, nmax);

        simdtrf::transform_g_inner(buffer, 60067, 55342, 15, 7, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 60067, 63, nmax);

        simdtrf::transform_g_inner(buffer, 60067, 56917, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 567 * nvalues + n * npairs, nvalues, buffer, 60067,
                                   63, nmax);

        simdtrf::transform_g_inner(buffer, 60067, 58492, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 1134 * nvalues + n * npairs, nvalues, buffer, 60067,
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
