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


#include "SimdThreeCenterElectronRepulsionGeom010RecPSH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
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
#include "SimdTransformH.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_psh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_psh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2515, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 99 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 2515, 2293, 189, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto beta = b_exps[j];

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 9, 6, {1, 2, 3, 4,
                                                        5, 6, 7}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 17, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 20, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 53, 0, 6, 10, 11,
                                                                       17, 20, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 62, 0, 6, 11, 12,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 71, 0, 6, 12, 13,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 80, 0, 6, 13, 14,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 89, 0, 6, 14, 15,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 98, 6, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 101, 6, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 104, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 107, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 110, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 113, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 116, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 119, 6, 12, 23,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 128, 6, 13, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 137, 6, 14, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 146, 6, 15, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 155, 6, 10, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 164, 6, 11, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 173, 6, 12, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 182, 6, 13, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 191, 6, 14, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 200, 6, 15, 50,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 209, 6, 17, 35,
                                                                       53, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 236, 6, 20, 38,
                                                                       62, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 263, 6, 23, 41,
                                                                       71, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 290, 6, 26, 44,
                                                                       80, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 317, 6, 29, 47,
                                                                       89, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 344, 6, 10, 11,
                                                                       104, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 350, 6, 11, 12,
                                                                       107, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 356, 6, 12, 13,
                                                                       110, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 362, 6, 13, 14,
                                                                       113, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 368, 6, 14, 15,
                                                                       116, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 374, 3, 6, 344,
                                                                       104, 350, 119, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 392, 3, 6, 350,
                                                                       107, 356, 128, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 410, 3, 6, 356,
                                                                       110, 362, 137, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 428, 3, 6, 362,
                                                                       113, 368, 146, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 446, 0, 6, 344,
                                                                       104, 350, 173, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 464, 0, 6, 350,
                                                                       107, 356, 182, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 482, 0, 6, 356,
                                                                       110, 362, 191, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 500, 0, 6, 362,
                                                                       113, 368, 200, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 518, 0, 3, 6, 374,
                                                                       119, 392, 446, 173, 464,
                                                                       53, 62, 263, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 572, 0, 3, 6, 392,
                                                                       128, 410, 464, 182, 482,
                                                                       62, 71, 290, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 626, 0, 3, 6, 410,
                                                                       137, 428, 482, 191, 500,
                                                                       71, 80, 317, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 680, 6, 98, 101,
                                                                       344, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 690, 6, 101, 104,
                                                                       350, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 700, 6, 104, 107,
                                                                       356, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 710, 6, 107, 110,
                                                                       362, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 720, 6, 110, 113,
                                                                       368, ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 730, 3, 6, 680,
                                                                       344, 690, 374, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 760, 3, 6, 690,
                                                                       350, 700, 392, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 790, 3, 6, 700,
                                                                       356, 710, 410, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 820, 3, 6, 710,
                                                                       362, 720, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 850, 0, 6, 680,
                                                                       344, 690, 155, 164, 446,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 880, 0, 6, 690,
                                                                       350, 700, 164, 173, 464,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 910, 0, 6, 700,
                                                                       356, 710, 173, 182, 482,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 940, 0, 6, 710,
                                                                       362, 720, 182, 191, 500,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 970, 0, 3, 6, 730,
                                                                       374, 760, 850, 446, 880,
                                                                       209, 236, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 1060, 0, 3, 6,
                                                                       760, 392, 790, 880, 464,
                                                                       910, 236, 263, 572, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 1150, 0, 3, 6,
                                                                       790, 410, 820, 910, 482,
                                                                       940, 263, 290, 626, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1240, 6, 344, 350,
                                                                       700, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1255, 6, 350, 356,
                                                                       710, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1270, 6, 356, 362,
                                                                       720, ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 1285, 3, 6, 1240,
                                                                       700, 1255, 374, 392, 790,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 1330, 3, 6, 1255,
                                                                       710, 1270, 392, 410, 820,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1375, 0, 6, 1240,
                                                                       700, 1255, 446, 464, 910,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1420, 0, 6, 1255,
                                                                       710, 1270, 464, 482, 940,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 1465, 0, 3, 6,
                                                                       1285, 790, 1330, 1375,
                                                                       910, 1420, 518, 572, 1150,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 1600, 6, 680, 690,
                                                                       1240, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 1621, 6, 690, 700,
                                                                       1255, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 1642, 6, 700, 710,
                                                                       1270, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 1663, 3, 6, 1600,
                                                                       1240, 1621, 730, 760,
                                                                       1285, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 1726, 3, 6, 1621,
                                                                       1255, 1642, 760, 790,
                                                                       1330, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 1789, 0, 6, 1600,
                                                                       1240, 1621, 850, 880,
                                                                       1375, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 1852, 0, 6, 1621,
                                                                       1255, 1642, 880, 910,
                                                                       1420, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 1915, 0, 3, 6,
                                                                       1663, 1285, 1726, 1789,
                                                                       1375, 1852, 970, 1060,
                                                                       1465, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_s_x(buffer, 2104, 1915, 3, 21, ncols, beta);

                    simdgeo::geom_s_y(buffer, 2167, 1915, 3, 21, ncols, beta);

                    simdgeo::geom_s_z(buffer, 2230, 1915, 3, 21, ncols, beta);

                    simdfunc::contract_primitives(buffer, 2293, 2104, 189, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 2482, 2293, 3, 1, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 2482, 11, nmax);

        simdtrf::transform_h_inner(buffer, 2482, 2356, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 33 * nvalues + n * npairs, nvalues, buffer, 2482, 11,
                                   nmax);

        simdtrf::transform_h_inner(buffer, 2482, 2419, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 66 * nvalues + n * npairs, nvalues, buffer, 2482, 11,
                                   nmax);
    }

    for (size_t m = 0; m < 99; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
