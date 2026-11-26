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


#include "SimdThreeCenterElectronRepulsionGeom100RecPFP.hpp"

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
#include "SimdGeometryF1.hpp"
#include "SimdGeometryG1.hpp"
#include "SimdGeometryP1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferGeom100XDD.hpp"
#include "SimdTransferGeom100XDP.hpp"
#include "SimdTransferGeom100XFP.hpp"
#include "SimdTransferGeom100XPD.hpp"
#include "SimdTransferGeom100XPF.hpp"
#include "SimdTransferGeom100XPP.hpp"
#include "SimdTransferGeom100YDD.hpp"
#include "SimdTransferGeom100YDP.hpp"
#include "SimdTransferGeom100YFP.hpp"
#include "SimdTransferGeom100YPD.hpp"
#include "SimdTransferGeom100YPF.hpp"
#include "SimdTransferGeom100YPP.hpp"
#include "SimdTransferGeom100ZDD.hpp"
#include "SimdTransferGeom100ZDP.hpp"
#include "SimdTransferGeom100ZFP.hpp"
#include "SimdTransferGeom100ZPD.hpp"
#include "SimdTransferGeom100ZPF.hpp"
#include "SimdTransferGeom100ZPP.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransferPP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_pfp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_pfp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2800, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 189 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 2800, 607, 681, dimensions);

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
                                                        5, 6}, ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 13, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 16, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 19, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 22, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 7, 8,
                                                                       13, 16, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 8, 9,
                                                                       16, 19, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 9, 10,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 10, 11,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 13, 16,
                                                                       28, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 16, 19,
                                                                       34, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 19, 22,
                                                                       40, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 28, 34,
                                                                       52, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 97, 0, 3, 34, 40,
                                                                       62, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 52, 62,
                                                                       82, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 133, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 136, 3, 7, 13,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 145, 3, 13, 28,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 163, 3, 28, 52,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 193, 3, 52, 82,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 238, 3, 82, 112,
                                                                       ncols, p, q);

                    simdgeo::geom_p_x(buffer, 301, 133, 145, 1, 3, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 310, 133, 145, 1, 3, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 319, 133, 145, 1, 3, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 328, 136, 163, 1, 3, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 346, 136, 163, 1, 3, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 364, 136, 163, 1, 3, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 382, 145, 193, 1, 3, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 412, 145, 193, 1, 3, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 442, 145, 193, 1, 3, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 472, 163, 238, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 517, 163, 238, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 562, 163, 238, 1, 3, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 607, 301, 9, ncols);

                    simdfunc::contract_primitives(buffer, 625, 310, 9, ncols);

                    simdfunc::contract_primitives(buffer, 643, 319, 9, ncols);

                    simdfunc::contract_primitives(buffer, 661, 136, 9, ncols);

                    simdfunc::contract_primitives(buffer, 679, 328, 18, ncols);

                    simdfunc::contract_primitives(buffer, 715, 346, 18, ncols);

                    simdfunc::contract_primitives(buffer, 751, 364, 18, ncols);

                    simdfunc::contract_primitives(buffer, 787, 145, 18, ncols);

                    simdfunc::contract_primitives(buffer, 823, 382, 30, ncols);

                    simdfunc::contract_primitives(buffer, 883, 412, 30, ncols);

                    simdfunc::contract_primitives(buffer, 943, 442, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1003, 163, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1063, 472, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1153, 517, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1243, 562, 45, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 616, 607, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 634, 625, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 652, 643, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 670, 661, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 697, 679, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 733, 715, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 769, 751, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 805, 787, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 853, 823, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 913, 883, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 973, 943, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1033, 1003, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1108, 1063, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1198, 1153, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1288, 1243, 15, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 1333, 616, 670, 697,
                                                       3, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 1360, 634, 670, 733,
                                                       3, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 1387, 652, 670, 769,
                                                       3, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 1414, 670, 805, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 1441, 697, 805, 853,
                                                       3, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 1495, 733, 805, 913,
                                                       3, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 1549, 769, 805, 973,
                                                       3, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 1603, 805, 1033, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 1657, 853, 1033,
                                                       1108, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 1747, 913, 1033,
                                                       1198, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 1837, 973, 1033,
                                                       1288, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 1927, 1333, 1414,
                                                       1441, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 1981, 1360, 1414,
                                                       1495, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 2035, 1387, 1414,
                                                       1549, 3, nmax);

        simdtrf::compute_hrr_pd_out_of_first(buffer, coordinates, 2089, 1414, 1603, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 2143, 1441, 1603,
                                                       1657, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 2251, 1495, 1603,
                                                       1747, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 2359, 1549, 1603,
                                                       1837, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pf_out_of_first(buffer, coordinates, 2467, 1927, 2089,
                                                       2143, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pf_out_of_first(buffer, coordinates, 2557, 1981, 2089,
                                                       2251, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pf_out_of_first(buffer, coordinates, 2647, 2035, 2089,
                                                       2359, 3, nmax);

        simdtrf::transform_f_inner(buffer, 2737, 2467, 3, 3, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 2737, 21, nmax);

        simdtrf::transform_f_inner(buffer, 2737, 2557, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 2737, 21,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 2737, 2647, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 126 * nvalues + n * npairs, nvalues, buffer, 2737,
                                   21, nmax);
    }

    for (size_t m = 0; m < 189; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
