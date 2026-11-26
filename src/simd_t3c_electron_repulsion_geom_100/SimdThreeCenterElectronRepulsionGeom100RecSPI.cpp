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


#include "SimdThreeCenterElectronRepulsionGeom100RecSPI.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdGeometryS1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
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
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_spi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_spi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 4035, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 117 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 4035, 3744, 252, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto alpha = a_exps[i];

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fa = -b_exps[j] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pa(buffer, coordinates, 0, nmax, fa);

                simdfunc::compute_pb(buffer, coordinates, 3, nmax, fb);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 6, n, nmax, fc);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 9, 6, 8, ncols,
                                                             fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 19, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 61, 0, 6, 10, 11,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 70, 0, 6, 11, 12,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 79, 0, 6, 12, 13,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 88, 0, 6, 13, 14,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 97, 0, 6, 14, 15,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 106, 0, 6, 15, 16,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 115, 0, 6, 16, 17,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 124, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 127, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 130, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 133, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 136, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 139, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 142, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 145, 6, 12, 25,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 154, 6, 13, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 163, 6, 14, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 172, 6, 15, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 181, 6, 16, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 190, 6, 17, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 199, 6, 12, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 208, 6, 13, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 217, 6, 14, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 226, 6, 15, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 235, 6, 16, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 244, 6, 17, 58,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 253, 6, 25, 43,
                                                                       79, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 280, 6, 28, 46,
                                                                       88, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 307, 6, 31, 49,
                                                                       97, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 334, 6, 34, 52,
                                                                       106, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 361, 6, 37, 55,
                                                                       115, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 388, 6, 10, 11,
                                                                       124, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 394, 6, 11, 12,
                                                                       127, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 400, 6, 12, 13,
                                                                       130, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 406, 6, 13, 14,
                                                                       133, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 412, 6, 14, 15,
                                                                       136, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 418, 6, 15, 16,
                                                                       139, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 424, 6, 16, 17,
                                                                       142, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 430, 3, 6, 388,
                                                                       124, 394, 145, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 448, 3, 6, 394,
                                                                       127, 400, 154, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 466, 3, 6, 400,
                                                                       130, 406, 163, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 484, 3, 6, 406,
                                                                       133, 412, 172, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 502, 3, 6, 412,
                                                                       136, 418, 181, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 520, 3, 6, 418,
                                                                       139, 424, 190, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 538, 0, 6, 388,
                                                                       124, 394, 199, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 556, 0, 6, 394,
                                                                       127, 400, 208, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 574, 0, 6, 400,
                                                                       130, 406, 217, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 592, 0, 6, 406,
                                                                       133, 412, 226, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 610, 0, 6, 412,
                                                                       136, 418, 235, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 628, 0, 6, 418,
                                                                       139, 424, 244, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 646, 0, 3, 6, 430,
                                                                       145, 448, 538, 199, 556,
                                                                       61, 70, 253, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 700, 0, 3, 6, 448,
                                                                       154, 466, 556, 208, 574,
                                                                       70, 79, 280, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 754, 0, 3, 6, 466,
                                                                       163, 484, 574, 217, 592,
                                                                       79, 88, 307, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 808, 0, 3, 6, 484,
                                                                       172, 502, 592, 226, 610,
                                                                       88, 97, 334, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 862, 0, 3, 6, 502,
                                                                       181, 520, 610, 235, 628,
                                                                       97, 106, 361, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 916, 6, 124, 127,
                                                                       400, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 926, 6, 127, 130,
                                                                       406, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 936, 6, 130, 133,
                                                                       412, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 946, 6, 133, 136,
                                                                       418, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 956, 6, 136, 139,
                                                                       424, ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 966, 3, 6, 916,
                                                                       400, 926, 466, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 996, 3, 6, 926,
                                                                       406, 936, 484, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1026, 3, 6, 936,
                                                                       412, 946, 502, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1056, 3, 6, 946,
                                                                       418, 956, 520, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1086, 0, 6, 916,
                                                                       400, 926, 199, 208, 574,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1116, 0, 6, 926,
                                                                       406, 936, 208, 217, 592,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1146, 0, 6, 936,
                                                                       412, 946, 217, 226, 610,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1176, 0, 6, 946,
                                                                       418, 956, 226, 235, 628,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 1206, 0, 3, 6,
                                                                       966, 466, 996, 1086, 574,
                                                                       1116, 253, 280, 754,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 6,
                                                                       996, 484, 1026, 1116, 592,
                                                                       1146, 280, 307, 808,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 1386, 0, 3, 6,
                                                                       1026, 502, 1056, 1146,
                                                                       610, 1176, 307, 334, 862,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1476, 6, 388, 394,
                                                                       916, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1491, 6, 394, 400,
                                                                       926, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1506, 6, 400, 406,
                                                                       936, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1521, 6, 406, 412,
                                                                       946, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1536, 6, 412, 418,
                                                                       956, ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 1551, 3, 6, 1476,
                                                                       916, 1491, 430, 448, 966,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 1596, 3, 6, 1491,
                                                                       926, 1506, 448, 466, 996,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 1641, 3, 6, 1506,
                                                                       936, 1521, 466, 484, 1026,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 1686, 3, 6, 1521,
                                                                       946, 1536, 484, 502, 1056,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1731, 0, 6, 1476,
                                                                       916, 1491, 538, 556, 1086,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1776, 0, 6, 1491,
                                                                       926, 1506, 556, 574, 1116,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1821, 0, 6, 1506,
                                                                       936, 1521, 574, 592, 1146,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1866, 0, 6, 1521,
                                                                       946, 1536, 592, 610, 1176,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 1911, 0, 3, 6,
                                                                       1551, 966, 1596, 1731,
                                                                       1086, 1776, 646, 700,
                                                                       1206, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 2046, 0, 3, 6,
                                                                       1596, 996, 1641, 1776,
                                                                       1116, 1821, 700, 754,
                                                                       1296, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 2181, 0, 3, 6,
                                                                       1641, 1026, 1686, 1821,
                                                                       1146, 1866, 754, 808,
                                                                       1386, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 2316, 6, 916, 926,
                                                                       1506, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 2337, 6, 926, 936,
                                                                       1521, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 2358, 6, 936, 946,
                                                                       1536, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 2379, 3, 6, 2316,
                                                                       1506, 2337, 966, 996,
                                                                       1641, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 2442, 3, 6, 2337,
                                                                       1521, 2358, 996, 1026,
                                                                       1686, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 2505, 0, 6, 2316,
                                                                       1506, 2337, 1086, 1116,
                                                                       1821, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 2568, 0, 6, 2337,
                                                                       1521, 2358, 1116, 1146,
                                                                       1866, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 2631, 0, 3, 6,
                                                                       2379, 1641, 2442, 2505,
                                                                       1821, 2568, 1206, 1296,
                                                                       2181, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 2820, 6, 1476,
                                                                       1491, 2316, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 2848, 6, 1491,
                                                                       1506, 2337, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 2876, 6, 1506,
                                                                       1521, 2358, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 2904, 3, 6, 2820,
                                                                       2316, 2848, 1551, 1596,
                                                                       2379, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 2988, 3, 6, 2848,
                                                                       2337, 2876, 1596, 1641,
                                                                       2442, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 3072, 0, 6, 2820,
                                                                       2316, 2848, 1731, 1776,
                                                                       2505, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 3156, 0, 6, 2848,
                                                                       2337, 2876, 1776, 1821,
                                                                       2568, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 3240, 0, 3, 6,
                                                                       2904, 2379, 2988, 3072,
                                                                       2505, 3156, 1911, 2046,
                                                                       2631, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_s_x(buffer, 3492, 3240, 1, 84, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 3576, 3240, 1, 84, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 3660, 3240, 1, 84, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 3744, 3492, 252, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 3996, 3744, 3, 1, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 3996, 13, nmax);

        simdtrf::transform_i_inner(buffer, 3996, 3828, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 39 * nvalues + n * npairs, nvalues, buffer, 3996, 13,
                                   nmax);

        simdtrf::transform_i_inner(buffer, 3996, 3912, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 78 * nvalues + n * npairs, nvalues, buffer, 3996, 13,
                                   nmax);
    }

    for (size_t m = 0; m < 117; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
