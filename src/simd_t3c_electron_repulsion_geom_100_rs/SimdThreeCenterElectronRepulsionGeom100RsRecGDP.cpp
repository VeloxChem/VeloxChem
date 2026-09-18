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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecGDP.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferGeom100XGD.hpp"
#include "SimdTransferGeom100XGP.hpp"
#include "SimdTransferGeom100XHP.hpp"
#include "SimdTransferGeom100YGD.hpp"
#include "SimdTransferGeom100YGP.hpp"
#include "SimdTransferGeom100YHP.hpp"
#include "SimdTransferGeom100ZGD.hpp"
#include "SimdTransferGeom100ZGP.hpp"
#include "SimdTransferGeom100ZHP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_gdp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_gdp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 9275, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 810 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 9275, 2480, 2652, dimensions);

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
                                                            4, 5, 6, 7, 8}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 15, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 24, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 27, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 60, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 63, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 66, 0, 3, 7, 8,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 8, 9,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 78, 0, 3, 9, 10,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 84, 0, 3, 10, 11,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 90, 0, 3, 11, 12,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 96, 0, 3, 12, 13,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 16, 17,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 17, 18,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 114, 0, 3, 18, 19,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 120, 0, 3, 19, 20,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 126, 0, 3, 20, 21,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 21, 22,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 24, 27,
                                                                       66, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 27, 30,
                                                                       72, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 30, 33,
                                                                       78, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 168, 0, 3, 33, 36,
                                                                       84, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 178, 0, 3, 36, 39,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 45, 48,
                                                                       102, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 198, 0, 3, 48, 51,
                                                                       108, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 208, 0, 3, 51, 54,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 54, 57,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 228, 0, 3, 57, 60,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 66, 72,
                                                                       138, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 253, 0, 3, 72, 78,
                                                                       148, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 78, 84,
                                                                       158, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 283, 0, 3, 84, 90,
                                                                       168, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 102,
                                                                       108, 188, 198, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 313, 0, 3, 108,
                                                                       114, 198, 208, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 114,
                                                                       120, 208, 218, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 343, 0, 3, 120,
                                                                       126, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 138,
                                                                       148, 238, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 379, 0, 3, 148,
                                                                       158, 253, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 400, 0, 3, 158,
                                                                       168, 268, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 421, 0, 3, 188,
                                                                       198, 298, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 442, 0, 3, 198,
                                                                       208, 313, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 463, 0, 3, 208,
                                                                       218, 328, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 484, 0, 3, 238,
                                                                       253, 358, 379, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 512, 0, 3, 253,
                                                                       268, 379, 400, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 540, 0, 3, 298,
                                                                       313, 421, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 568, 0, 3, 313,
                                                                       328, 442, 463, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 596, 0, 3, 358,
                                                                       379, 484, 512, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 632, 0, 3, 421,
                                                                       442, 540, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 668, 3, 66, 138,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 698, 3, 102, 188,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 728, 3, 138, 238,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 773, 3, 188, 298,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 818, 3, 238, 358,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 881, 3, 298, 421,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 944, 3, 358, 484,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 1028, 3, 421, 540,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 1112, 3, 484, 596,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 1220, 3, 540, 632,
                                                                       ncols, p, q);

                    simdgeo::geom_g_x(buffer, 1328, 668, 818, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 1373, 668, 818, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 1418, 668, 818, 1, 3, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 1463, 698, 881, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 1508, 698, 881, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 1553, 698, 881, 1, 3, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 1598, 728, 944, 1, 3, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 1661, 728, 944, 1, 3, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 1724, 728, 944, 1, 3, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 1787, 773, 1028, 1, 3, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 1850, 773, 1028, 1, 3, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 1913, 773, 1028, 1, 3, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 1976, 818, 1112, 1, 3, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 2060, 818, 1112, 1, 3, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 2144, 818, 1112, 1, 3, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 2228, 881, 1220, 1, 3, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 2312, 881, 1220, 1, 3, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 2396, 881, 1220, 1, 3, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 2480, 1328, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2570, 1373, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2660, 1418, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2750, 728, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2840, 1463, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2930, 1508, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3020, 1553, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3110, 773, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3200, 1598, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3326, 1661, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3452, 1724, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3578, 818, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3704, 1787, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3830, 1850, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3956, 1913, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4082, 881, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4208, 1976, 84, ncols);

                    simdfunc::contract_primitives(buffer, 4376, 2060, 84, ncols);

                    simdfunc::contract_primitives(buffer, 4544, 2144, 84, ncols);

                    simdfunc::contract_primitives(buffer, 4712, 2228, 84, ncols);

                    simdfunc::contract_primitives(buffer, 4880, 2312, 84, ncols);

                    simdfunc::contract_primitives(buffer, 5048, 2396, 84, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 2525, 2480, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2615, 2570, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2705, 2660, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2795, 2750, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2885, 2840, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2975, 2930, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3065, 3020, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3155, 3110, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3263, 3200, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3389, 3326, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3515, 3452, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3641, 3578, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3767, 3704, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3893, 3830, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4019, 3956, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4145, 4082, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4292, 4208, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4460, 4376, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4628, 4544, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4796, 4712, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4964, 4880, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5132, 5048, 28, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 5216, 2525, 2795,
                                                       3263, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 5351, 2615, 2795,
                                                       3389, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 5486, 2705, 2795,
                                                       3515, 3, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 5621, 2795, 3641, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 5756, 2885, 3155,
                                                       3767, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 5891, 2975, 3155,
                                                       3893, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 6026, 3065, 3155,
                                                       4019, 3, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 6161, 3155, 4145, 3, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 6296, 3263, 3641,
                                                       4292, 3, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 6485, 3389, 3641,
                                                       4460, 3, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 6674, 3515, 3641,
                                                       4628, 3, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 6863, 3767, 4145,
                                                       4796, 3, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 7052, 3893, 4145,
                                                       4964, 3, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 7241, 4019, 4145,
                                                       5132, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 7430, 5216, 5621,
                                                       6296, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 7700, 5351, 5621,
                                                       6485, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 7970, 5486, 5621,
                                                       6674, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 8240, 5756, 6161,
                                                       6863, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 8510, 5891, 6161,
                                                       7052, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 8780, 6026, 6161,
                                                       7241, 3, nmax);

        simdtrf::transform_d_inner(buffer, 9050, 8240, 15, 3, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 9050, 15, nmax);

        simdtrf::transform_d_inner(buffer, 9050, 8510, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 135 * nvalues + n * npairs, nvalues, buffer, 9050,
                                   15, nmax);

        simdtrf::transform_d_inner(buffer, 9050, 8780, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 270 * nvalues + n * npairs, nvalues, buffer, 9050,
                                   15, nmax);

        simdtrf::transform_d_inner(buffer, 9050, 7430, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 405 * nvalues + n * npairs, nvalues, buffer, 9050,
                                   15, nmax);

        simdtrf::transform_d_inner(buffer, 9050, 7700, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 540 * nvalues + n * npairs, nvalues, buffer, 9050,
                                   15, nmax);

        simdtrf::transform_d_inner(buffer, 9050, 7970, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 675 * nvalues + n * npairs, nvalues, buffer, 9050,
                                   15, nmax);
    }

    for (size_t m = 0; m < 810; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
