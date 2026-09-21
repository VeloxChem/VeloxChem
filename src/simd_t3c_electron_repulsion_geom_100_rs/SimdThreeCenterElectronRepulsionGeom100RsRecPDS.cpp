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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecPDS.hpp"

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
#include "SimdGeometryP1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdTransferGeom100XDP.hpp"
#include "SimdTransferGeom100XPD.hpp"
#include "SimdTransferGeom100XPP.hpp"
#include "SimdTransferGeom100YDP.hpp"
#include "SimdTransferGeom100YPD.hpp"
#include "SimdTransferGeom100YPP.hpp"
#include "SimdTransferGeom100ZDP.hpp"
#include "SimdTransferGeom100ZPD.hpp"
#include "SimdTransferGeom100ZPP.hpp"
#include "SimdTransferPP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformP.hpp"
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_pds_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_pds_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 829, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 90 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 829, 262, 254, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 4,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 12, 3, 4,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 18, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 21, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 24, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 27, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 7, 8,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 8, 9,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 9, 10,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 60, 0, 3, 13, 14,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 66, 0, 3, 14, 15,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 15, 16,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 78, 0, 3, 18, 21,
                                                                       42, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 88, 0, 3, 21, 24,
                                                                       48, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 30, 33,
                                                                       60, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 33, 36,
                                                                       66, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 42, 48,
                                                                       78, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 133, 0, 3, 60, 66,
                                                                       98, 108, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_p_x(buffer, 148, 7, 42, 1, 1, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 151, 7, 42, 1, 1, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 154, 7, 42, 1, 1, ncols, alpha);

                    simdgeo::geom_p_x(buffer, 157, 13, 60, 1, 1, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 160, 13, 60, 1, 1, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 163, 13, 60, 1, 1, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 166, 18, 78, 1, 1, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 172, 18, 78, 1, 1, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 178, 18, 78, 1, 1, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 184, 30, 98, 1, 1, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 190, 30, 98, 1, 1, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 196, 30, 98, 1, 1, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 202, 42, 118, 1, 1, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 212, 42, 118, 1, 1, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 222, 42, 118, 1, 1, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 232, 60, 133, 1, 1, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 242, 60, 133, 1, 1, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 252, 60, 133, 1, 1, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 262, 148, 3, ncols);

                    simdfunc::contract_primitives(buffer, 268, 151, 3, ncols);

                    simdfunc::contract_primitives(buffer, 274, 154, 3, ncols);

                    simdfunc::contract_primitives(buffer, 280, 18, 3, ncols);

                    simdfunc::contract_primitives(buffer, 286, 157, 3, ncols);

                    simdfunc::contract_primitives(buffer, 292, 160, 3, ncols);

                    simdfunc::contract_primitives(buffer, 298, 163, 3, ncols);

                    simdfunc::contract_primitives(buffer, 304, 30, 3, ncols);

                    simdfunc::contract_primitives(buffer, 310, 166, 6, ncols);

                    simdfunc::contract_primitives(buffer, 322, 172, 6, ncols);

                    simdfunc::contract_primitives(buffer, 334, 178, 6, ncols);

                    simdfunc::contract_primitives(buffer, 346, 42, 6, ncols);

                    simdfunc::contract_primitives(buffer, 358, 184, 6, ncols);

                    simdfunc::contract_primitives(buffer, 370, 190, 6, ncols);

                    simdfunc::contract_primitives(buffer, 382, 196, 6, ncols);

                    simdfunc::contract_primitives(buffer, 394, 60, 6, ncols);

                    simdfunc::contract_primitives(buffer, 406, 202, 10, ncols);

                    simdfunc::contract_primitives(buffer, 426, 212, 10, ncols);

                    simdfunc::contract_primitives(buffer, 446, 222, 10, ncols);

                    simdfunc::contract_primitives(buffer, 466, 232, 10, ncols);

                    simdfunc::contract_primitives(buffer, 486, 242, 10, ncols);

                    simdfunc::contract_primitives(buffer, 506, 252, 10, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 265, 262, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 271, 268, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 277, 274, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 283, 280, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 289, 286, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 295, 292, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 301, 298, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 307, 304, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 316, 310, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 328, 322, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 340, 334, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 352, 346, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 364, 358, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 376, 370, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 388, 382, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 400, 394, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 416, 406, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 436, 426, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 456, 446, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 476, 466, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 496, 486, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 516, 506, 10, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 526, 265, 283, 316,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 535, 271, 283, 328,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 544, 277, 283, 340,
                                                       1, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 553, 283, 352, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 562, 289, 307, 364,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 571, 295, 307, 376,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 580, 301, 307, 388,
                                                       1, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 589, 307, 400, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 598, 316, 352, 416,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 616, 328, 352, 436,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 634, 340, 352, 456,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 652, 364, 400, 476,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 670, 376, 400, 496,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 688, 388, 400, 516,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 706, 526, 553, 598,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 724, 535, 553, 616,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 742, 544, 553, 634,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 760, 562, 589, 652,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 778, 571, 589, 670,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 796, 580, 589, 688,
                                                       1, nmax);

        simdtrf::transform_d_inner(buffer, 814, 760, 3, 1, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 814, 5, nmax);

        simdtrf::transform_d_inner(buffer, 814, 778, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 15 * nvalues + n * npairs, nvalues, buffer, 814, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 814, 796, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 30 * nvalues + n * npairs, nvalues, buffer, 814, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 814, 706, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 45 * nvalues + n * npairs, nvalues, buffer, 814, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 814, 724, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 60 * nvalues + n * npairs, nvalues, buffer, 814, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 814, 742, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 75 * nvalues + n * npairs, nvalues, buffer, 814, 5,
                                   nmax);
    }

    for (size_t m = 0; m < 90; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
