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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecPPF.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdGeometryD1.hpp"
#include "SimdGeometryP1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferGeom100XPP.hpp"
#include "SimdTransferGeom100YPP.hpp"
#include "SimdTransferGeom100ZPP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_ppf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_ppf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 3729, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 378 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 3729, 2268, 978, dimensions);

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
                                                            4, 5, 6}, ncols, fj, i * nprim_b + j,
                                                            fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 13, 3, {1, 2, 3, 4,
                                                        5, 6}, ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 20, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 23, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 7, 8,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 8, 9,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 9, 10,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 10, 11,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 14, 15,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 15, 16,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 16, 17,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 17, 18,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 20, 23,
                                                                       50, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 23, 26,
                                                                       56, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 26, 29,
                                                                       62, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 35, 38,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 38, 41,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 41, 44,
                                                                       86, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 158, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 161, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 164, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 167, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 170, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 173, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 176, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 179, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 182, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 185, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 188, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 191, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 194, 3, 7, 20,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 203, 3, 8, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 212, 3, 9, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 221, 3, 10, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 230, 3, 11, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 239, 3, 14, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 248, 3, 15, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 257, 3, 16, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 266, 3, 17, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 275, 3, 18, 47,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 284, 3, 20, 50,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 302, 3, 23, 56,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 320, 3, 26, 62,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 338, 3, 29, 68,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 356, 3, 35, 74,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 374, 3, 38, 80,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 392, 3, 41, 86,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 410, 3, 44, 92,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 428, 3, 50, 98,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 458, 3, 56, 108,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 488, 3, 62, 118,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 518, 3, 74, 128,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 548, 3, 80, 138,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 578, 3, 86, 148,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 608, 3, 7, 8, 164,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 614, 3, 8, 9, 167,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 620, 3, 9, 10,
                                                                       170, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 626, 3, 10, 11,
                                                                       173, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 632, 3, 14, 15,
                                                                       182, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 638, 3, 15, 16,
                                                                       185, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 644, 3, 16, 17,
                                                                       188, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 650, 3, 17, 18,
                                                                       191, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 656, 0, 3, 608,
                                                                       164, 614, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 674, 0, 3, 614,
                                                                       167, 620, 221, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 692, 0, 3, 620,
                                                                       170, 626, 230, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 710, 0, 3, 632,
                                                                       182, 638, 257, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 728, 0, 3, 638,
                                                                       185, 644, 266, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 746, 0, 3, 644,
                                                                       188, 650, 275, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 764, 0, 3, 656,
                                                                       212, 674, 50, 56, 320,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 800, 0, 3, 674,
                                                                       221, 692, 56, 62, 338,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 836, 0, 3, 710,
                                                                       257, 728, 74, 80, 392,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 872, 0, 3, 728,
                                                                       266, 746, 80, 86, 410,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 908, 0, 3, 764,
                                                                       320, 800, 98, 108, 488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 968, 0, 3, 836,
                                                                       392, 872, 128, 138, 578,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1028, 3, 158, 161,
                                                                       608, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1038, 3, 161, 164,
                                                                       614, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1048, 3, 164, 167,
                                                                       620, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1058, 3, 167, 170,
                                                                       626, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1068, 3, 176, 179,
                                                                       632, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1078, 3, 179, 182,
                                                                       638, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1088, 3, 182, 185,
                                                                       644, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1098, 3, 185, 188,
                                                                       650, ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1108, 0, 3, 1028,
                                                                       608, 1038, 194, 203, 656,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1138, 0, 3, 1038,
                                                                       614, 1048, 203, 212, 674,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1168, 0, 3, 1048,
                                                                       620, 1058, 212, 221, 692,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1198, 0, 3, 1068,
                                                                       632, 1078, 239, 248, 710,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1228, 0, 3, 1078,
                                                                       638, 1088, 248, 257, 728,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1258, 0, 3, 1088,
                                                                       644, 1098, 257, 266, 746,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1288, 0, 3, 1108,
                                                                       656, 1138, 284, 302, 764,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1348, 0, 3, 1138,
                                                                       674, 1168, 302, 320, 800,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 1198,
                                                                       710, 1228, 356, 374, 836,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1468, 0, 3, 1228,
                                                                       728, 1258, 374, 392, 872,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 1528, 0, 3, 1288,
                                                                       764, 1348, 428, 458, 908,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 1628, 0, 3, 1408,
                                                                       836, 1468, 518, 548, 968,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_p_x(buffer, 1728, 1028, 1288, 1, 10, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 1758, 1028, 1288, 1, 10, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 1788, 1028, 1288, 1, 10, ncols, alpha);

                    simdgeo::geom_p_x(buffer, 1818, 1068, 1408, 1, 10, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 1848, 1068, 1408, 1, 10, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 1878, 1068, 1408, 1, 10, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 1908, 1108, 1528, 1, 10, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 1968, 1108, 1528, 1, 10, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 2028, 1108, 1528, 1, 10, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 2088, 1198, 1628, 1, 10, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 2148, 1198, 1628, 1, 10, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 2208, 1198, 1628, 1, 10, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 2268, 1728, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2319, 1758, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2370, 1788, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2421, 1108, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2472, 1818, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2523, 1848, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2574, 1878, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2625, 1198, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2676, 1908, 60, ncols);

                    simdfunc::contract_primitives(buffer, 2778, 1968, 60, ncols);

                    simdfunc::contract_primitives(buffer, 2880, 2028, 60, ncols);

                    simdfunc::contract_primitives(buffer, 2982, 2088, 60, ncols);

                    simdfunc::contract_primitives(buffer, 3084, 2148, 60, ncols);

                    simdfunc::contract_primitives(buffer, 3186, 2208, 60, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 2298, 2268, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2349, 2319, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2400, 2370, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2451, 2421, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2502, 2472, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2553, 2523, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2604, 2574, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2655, 2625, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2736, 2676, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2838, 2778, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2940, 2880, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 3042, 2982, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 3144, 3084, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 3246, 3186, 6, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 3288, 2298, 2451,
                                                       2736, 7, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 3351, 2349, 2451,
                                                       2838, 7, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 3414, 2400, 2451,
                                                       2940, 7, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 3477, 2502, 2655,
                                                       3042, 7, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 3540, 2553, 2655,
                                                       3144, 7, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 3603, 2604, 2655,
                                                       3246, 7, nmax);

        simdtrf::transform_p_inner(buffer, 3666, 3477, 3, 7, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 3666, 21, nmax);

        simdtrf::transform_p_inner(buffer, 3666, 3540, 3, 7, nmax);

        simdtrf::transform_p_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 3666, 21,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 3666, 3603, 3, 7, nmax);

        simdtrf::transform_p_outer(values + 126 * nvalues + n * npairs, nvalues, buffer, 3666,
                                   21, nmax);

        simdtrf::transform_p_inner(buffer, 3666, 3288, 3, 7, nmax);

        simdtrf::transform_p_outer(values + 189 * nvalues + n * npairs, nvalues, buffer, 3666,
                                   21, nmax);

        simdtrf::transform_p_inner(buffer, 3666, 3351, 3, 7, nmax);

        simdtrf::transform_p_outer(values + 252 * nvalues + n * npairs, nvalues, buffer, 3666,
                                   21, nmax);

        simdtrf::transform_p_inner(buffer, 3666, 3414, 3, 7, nmax);

        simdtrf::transform_p_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 3666,
                                   21, nmax);
    }

    for (size_t m = 0; m < 378; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
