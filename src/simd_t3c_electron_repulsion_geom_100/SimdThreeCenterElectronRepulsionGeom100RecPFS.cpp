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


#include "SimdThreeCenterElectronRepulsionGeom100RecPFS.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
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
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_pfs_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_pfs_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 966, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 63 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 966, 235, 227, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 5, ncols,
                                                             fj, i * nprim_b + j, fq);

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

                    simdgeo::geom_p_x(buffer, 133, 7, 28, 1, 1, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 136, 7, 28, 1, 1, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 139, 7, 28, 1, 1, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 142, 13, 52, 1, 1, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 148, 13, 52, 1, 1, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 154, 13, 52, 1, 1, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 160, 28, 82, 1, 1, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 170, 28, 82, 1, 1, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 180, 28, 82, 1, 1, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 190, 52, 112, 1, 1, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 205, 52, 112, 1, 1, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 220, 52, 112, 1, 1, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 235, 133, 3, ncols);

                    simdfunc::contract_primitives(buffer, 241, 136, 3, ncols);

                    simdfunc::contract_primitives(buffer, 247, 139, 3, ncols);

                    simdfunc::contract_primitives(buffer, 253, 13, 3, ncols);

                    simdfunc::contract_primitives(buffer, 259, 142, 6, ncols);

                    simdfunc::contract_primitives(buffer, 271, 148, 6, ncols);

                    simdfunc::contract_primitives(buffer, 283, 154, 6, ncols);

                    simdfunc::contract_primitives(buffer, 295, 28, 6, ncols);

                    simdfunc::contract_primitives(buffer, 307, 160, 10, ncols);

                    simdfunc::contract_primitives(buffer, 327, 170, 10, ncols);

                    simdfunc::contract_primitives(buffer, 347, 180, 10, ncols);

                    simdfunc::contract_primitives(buffer, 367, 52, 10, ncols);

                    simdfunc::contract_primitives(buffer, 387, 190, 15, ncols);

                    simdfunc::contract_primitives(buffer, 417, 205, 15, ncols);

                    simdfunc::contract_primitives(buffer, 447, 220, 15, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 238, 235, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 244, 241, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 250, 247, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 256, 253, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 265, 259, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 277, 271, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 289, 283, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 301, 295, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 317, 307, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 337, 327, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 357, 347, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 377, 367, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 402, 387, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 432, 417, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 462, 447, 15, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 477, 238, 256, 265,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 486, 244, 256, 277,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 495, 250, 256, 289,
                                                       1, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 504, 256, 301, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 513, 265, 301, 317,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 531, 277, 301, 337,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 549, 289, 301, 357,
                                                       1, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 567, 301, 377, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 585, 317, 377, 402,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 615, 337, 377, 432,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 645, 357, 377, 462,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 675, 477, 504, 513,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 693, 486, 504, 531,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 711, 495, 504, 549,
                                                       1, nmax);

        simdtrf::compute_hrr_pd_out_of_first(buffer, coordinates, 729, 504, 567, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 747, 513, 567, 585,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 783, 531, 567, 615,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 819, 549, 567, 645,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100x_pf_out_of_first(buffer, coordinates, 855, 675, 729, 747,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_pf_out_of_first(buffer, coordinates, 885, 693, 729, 783,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_pf_out_of_first(buffer, coordinates, 915, 711, 729, 819,
                                                       1, nmax);

        simdtrf::transform_f_inner(buffer, 945, 855, 3, 1, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 945, 7, nmax);

        simdtrf::transform_f_inner(buffer, 945, 885, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 21 * nvalues + n * npairs, nvalues, buffer, 945, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 945, 915, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 42 * nvalues + n * npairs, nvalues, buffer, 945, 7,
                                   nmax);
    }

    for (size_t m = 0; m < 63; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
