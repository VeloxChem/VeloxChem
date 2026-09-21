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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecPGS.hpp"

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
#include "SimdGeometryH1.hpp"
#include "SimdGeometryP1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdTransferDD.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGeom100XDD.hpp"
#include "SimdTransferGeom100XDF.hpp"
#include "SimdTransferGeom100XDP.hpp"
#include "SimdTransferGeom100XFD.hpp"
#include "SimdTransferGeom100XFP.hpp"
#include "SimdTransferGeom100XGP.hpp"
#include "SimdTransferGeom100XPD.hpp"
#include "SimdTransferGeom100XPF.hpp"
#include "SimdTransferGeom100XPG.hpp"
#include "SimdTransferGeom100XPP.hpp"
#include "SimdTransferGeom100YDD.hpp"
#include "SimdTransferGeom100YDF.hpp"
#include "SimdTransferGeom100YDP.hpp"
#include "SimdTransferGeom100YFD.hpp"
#include "SimdTransferGeom100YFP.hpp"
#include "SimdTransferGeom100YGP.hpp"
#include "SimdTransferGeom100YPD.hpp"
#include "SimdTransferGeom100YPF.hpp"
#include "SimdTransferGeom100YPG.hpp"
#include "SimdTransferGeom100YPP.hpp"
#include "SimdTransferGeom100ZDD.hpp"
#include "SimdTransferGeom100ZDF.hpp"
#include "SimdTransferGeom100ZDP.hpp"
#include "SimdTransferGeom100ZFD.hpp"
#include "SimdTransferGeom100ZFP.hpp"
#include "SimdTransferGeom100ZGP.hpp"
#include "SimdTransferGeom100ZPD.hpp"
#include "SimdTransferGeom100ZPF.hpp"
#include "SimdTransferGeom100ZPG.hpp"
#include "SimdTransferGeom100ZPP.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransferPP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_pgs_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_pgs_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 3969, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 162 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 3969, 758, 775, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 6,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 14, 3, 6,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 22, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 7, 8,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 8, 9,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 9, 10,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 76, 0, 3, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 88, 0, 3, 15, 16,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 94, 0, 3, 16, 17,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 100, 0, 3, 17, 18,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 106, 0, 3, 18, 19,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 19, 20,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 22, 25,
                                                                       58, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 25, 28,
                                                                       64, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 28, 31,
                                                                       70, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 31, 34,
                                                                       76, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 40, 43,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 168, 0, 3, 43, 46,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 178, 0, 3, 46, 49,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 49, 52,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 198, 0, 3, 58, 64,
                                                                       118, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 64, 70,
                                                                       128, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 228, 0, 3, 70, 76,
                                                                       138, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 243, 0, 3, 88, 94,
                                                                       158, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 94,
                                                                       100, 168, 178, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 100,
                                                                       106, 178, 188, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 118,
                                                                       128, 198, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 309, 0, 3, 128,
                                                                       138, 213, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 330, 0, 3, 158,
                                                                       168, 243, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 351, 0, 3, 168,
                                                                       178, 258, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 372, 0, 3, 198,
                                                                       213, 288, 309, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 400, 0, 3, 243,
                                                                       258, 330, 351, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_p_x(buffer, 428, 7, 58, 1, 1, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 431, 7, 58, 1, 1, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 434, 7, 58, 1, 1, ncols, alpha);

                    simdgeo::geom_p_x(buffer, 437, 15, 88, 1, 1, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 440, 15, 88, 1, 1, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 443, 15, 88, 1, 1, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 446, 22, 118, 1, 1, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 452, 22, 118, 1, 1, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 458, 22, 118, 1, 1, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 464, 40, 158, 1, 1, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 470, 40, 158, 1, 1, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 476, 40, 158, 1, 1, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 482, 58, 198, 1, 1, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 492, 58, 198, 1, 1, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 502, 58, 198, 1, 1, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 512, 88, 243, 1, 1, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 522, 88, 243, 1, 1, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 532, 88, 243, 1, 1, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 542, 118, 288, 1, 1, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 557, 118, 288, 1, 1, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 572, 118, 288, 1, 1, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 587, 158, 330, 1, 1, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 602, 158, 330, 1, 1, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 617, 158, 330, 1, 1, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 632, 198, 372, 1, 1, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 653, 198, 372, 1, 1, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 674, 198, 372, 1, 1, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 695, 243, 400, 1, 1, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 716, 243, 400, 1, 1, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 737, 243, 400, 1, 1, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 758, 428, 3, ncols);

                    simdfunc::contract_primitives(buffer, 764, 431, 3, ncols);

                    simdfunc::contract_primitives(buffer, 770, 434, 3, ncols);

                    simdfunc::contract_primitives(buffer, 776, 22, 3, ncols);

                    simdfunc::contract_primitives(buffer, 782, 437, 3, ncols);

                    simdfunc::contract_primitives(buffer, 788, 440, 3, ncols);

                    simdfunc::contract_primitives(buffer, 794, 443, 3, ncols);

                    simdfunc::contract_primitives(buffer, 800, 40, 3, ncols);

                    simdfunc::contract_primitives(buffer, 806, 446, 6, ncols);

                    simdfunc::contract_primitives(buffer, 818, 452, 6, ncols);

                    simdfunc::contract_primitives(buffer, 830, 458, 6, ncols);

                    simdfunc::contract_primitives(buffer, 842, 58, 6, ncols);

                    simdfunc::contract_primitives(buffer, 854, 464, 6, ncols);

                    simdfunc::contract_primitives(buffer, 866, 470, 6, ncols);

                    simdfunc::contract_primitives(buffer, 878, 476, 6, ncols);

                    simdfunc::contract_primitives(buffer, 890, 88, 6, ncols);

                    simdfunc::contract_primitives(buffer, 902, 482, 10, ncols);

                    simdfunc::contract_primitives(buffer, 922, 492, 10, ncols);

                    simdfunc::contract_primitives(buffer, 942, 502, 10, ncols);

                    simdfunc::contract_primitives(buffer, 962, 118, 10, ncols);

                    simdfunc::contract_primitives(buffer, 982, 512, 10, ncols);

                    simdfunc::contract_primitives(buffer, 1002, 522, 10, ncols);

                    simdfunc::contract_primitives(buffer, 1022, 532, 10, ncols);

                    simdfunc::contract_primitives(buffer, 1042, 158, 10, ncols);

                    simdfunc::contract_primitives(buffer, 1062, 542, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1092, 557, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1122, 572, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1152, 198, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1182, 587, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1212, 602, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1242, 617, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1272, 243, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1302, 632, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1344, 653, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1386, 674, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1428, 695, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1470, 716, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1512, 737, 21, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 761, 758, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 767, 764, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 773, 770, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 779, 776, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 785, 782, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 791, 788, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 797, 794, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 803, 800, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 812, 806, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 824, 818, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 836, 830, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 848, 842, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 860, 854, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 872, 866, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 884, 878, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 896, 890, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 912, 902, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 932, 922, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 952, 942, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 972, 962, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 992, 982, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1012, 1002, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1032, 1022, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1052, 1042, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1077, 1062, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1107, 1092, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1137, 1122, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1167, 1152, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1197, 1182, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1227, 1212, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1257, 1242, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1287, 1272, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1323, 1302, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1365, 1344, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1407, 1386, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1449, 1428, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1491, 1470, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1533, 1512, 21, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 1554, 761, 779, 812,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 1563, 767, 779, 824,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 1572, 773, 779, 836,
                                                       1, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 1581, 779, 848, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 1590, 785, 803, 860,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 1599, 791, 803, 872,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 1608, 797, 803, 884,
                                                       1, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 1617, 803, 896, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 1626, 812, 848, 912,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 1644, 824, 848, 932,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 1662, 836, 848, 952,
                                                       1, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 1680, 848, 972, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 1698, 860, 896, 992,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 1716, 872, 896, 1012,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 1734, 884, 896, 1032,
                                                       1, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 1752, 896, 1052, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 1770, 912, 972, 1077,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 1800, 932, 972, 1107,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 1830, 952, 972, 1137,
                                                       1, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 1860, 972, 1167, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 1890, 992, 1052,
                                                       1197, 1, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 1920, 1012, 1052,
                                                       1227, 1, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 1950, 1032, 1052,
                                                       1257, 1, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 1980, 1052, 1287, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 2010, 1077, 1167,
                                                       1323, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 2055, 1107, 1167,
                                                       1365, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 2100, 1137, 1167,
                                                       1407, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 2145, 1197, 1287,
                                                       1449, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 2190, 1227, 1287,
                                                       1491, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 2235, 1257, 1287,
                                                       1533, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 2280, 1554, 1581,
                                                       1626, 1, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 2298, 1563, 1581,
                                                       1644, 1, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 2316, 1572, 1581,
                                                       1662, 1, nmax);

        simdtrf::compute_hrr_pd_out_of_first(buffer, coordinates, 2334, 1581, 1680, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 2352, 1590, 1617,
                                                       1698, 1, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 2370, 1599, 1617,
                                                       1716, 1, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 2388, 1608, 1617,
                                                       1734, 1, nmax);

        simdtrf::compute_hrr_pd_out_of_first(buffer, coordinates, 2406, 1617, 1752, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 2424, 1626, 1680,
                                                       1770, 1, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 2460, 1644, 1680,
                                                       1800, 1, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 2496, 1662, 1680,
                                                       1830, 1, nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 2532, 1680, 1860, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 2568, 1698, 1752,
                                                       1890, 1, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 2604, 1716, 1752,
                                                       1920, 1, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 2640, 1734, 1752,
                                                       1950, 1, nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 2676, 1752, 1980, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 2712, 1770, 1860,
                                                       2010, 1, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 2772, 1800, 1860,
                                                       2055, 1, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 2832, 1830, 1860,
                                                       2100, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 2892, 1890, 1980,
                                                       2145, 1, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 2952, 1920, 1980,
                                                       2190, 1, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 3012, 1950, 1980,
                                                       2235, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pf_out_of_first(buffer, coordinates, 3072, 2280, 2334,
                                                       2424, 1, nmax);

        simdtrf::compute_hrr_geom_100y_pf_out_of_first(buffer, coordinates, 3102, 2298, 2334,
                                                       2460, 1, nmax);

        simdtrf::compute_hrr_geom_100z_pf_out_of_first(buffer, coordinates, 3132, 2316, 2334,
                                                       2496, 1, nmax);

        simdtrf::compute_hrr_pf_out_of_first(buffer, coordinates, 3162, 2334, 2532, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pf_out_of_first(buffer, coordinates, 3192, 2352, 2406,
                                                       2568, 1, nmax);

        simdtrf::compute_hrr_geom_100y_pf_out_of_first(buffer, coordinates, 3222, 2370, 2406,
                                                       2604, 1, nmax);

        simdtrf::compute_hrr_geom_100z_pf_out_of_first(buffer, coordinates, 3252, 2388, 2406,
                                                       2640, 1, nmax);

        simdtrf::compute_hrr_pf_out_of_first(buffer, coordinates, 3282, 2406, 2676, 1, nmax);

        simdtrf::compute_hrr_geom_100x_df_out_of_first(buffer, coordinates, 3312, 2424, 2532,
                                                       2712, 1, nmax);

        simdtrf::compute_hrr_geom_100y_df_out_of_first(buffer, coordinates, 3372, 2460, 2532,
                                                       2772, 1, nmax);

        simdtrf::compute_hrr_geom_100z_df_out_of_first(buffer, coordinates, 3432, 2496, 2532,
                                                       2832, 1, nmax);

        simdtrf::compute_hrr_geom_100x_df_out_of_first(buffer, coordinates, 3492, 2568, 2676,
                                                       2892, 1, nmax);

        simdtrf::compute_hrr_geom_100y_df_out_of_first(buffer, coordinates, 3552, 2604, 2676,
                                                       2952, 1, nmax);

        simdtrf::compute_hrr_geom_100z_df_out_of_first(buffer, coordinates, 3612, 2640, 2676,
                                                       3012, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pg_out_of_first(buffer, coordinates, 3672, 3072, 3162,
                                                       3312, 1, nmax);

        simdtrf::compute_hrr_geom_100y_pg_out_of_first(buffer, coordinates, 3717, 3102, 3162,
                                                       3372, 1, nmax);

        simdtrf::compute_hrr_geom_100z_pg_out_of_first(buffer, coordinates, 3762, 3132, 3162,
                                                       3432, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pg_out_of_first(buffer, coordinates, 3807, 3192, 3282,
                                                       3492, 1, nmax);

        simdtrf::compute_hrr_geom_100y_pg_out_of_first(buffer, coordinates, 3852, 3222, 3282,
                                                       3552, 1, nmax);

        simdtrf::compute_hrr_geom_100z_pg_out_of_first(buffer, coordinates, 3897, 3252, 3282,
                                                       3612, 1, nmax);

        simdtrf::transform_g_inner(buffer, 3942, 3807, 3, 1, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 3942, 9, nmax);

        simdtrf::transform_g_inner(buffer, 3942, 3852, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 27 * nvalues + n * npairs, nvalues, buffer, 3942, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 3942, 3897, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 54 * nvalues + n * npairs, nvalues, buffer, 3942, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 3942, 3672, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 81 * nvalues + n * npairs, nvalues, buffer, 3942, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 3942, 3717, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 108 * nvalues + n * npairs, nvalues, buffer, 3942, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 3942, 3762, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 135 * nvalues + n * npairs, nvalues, buffer, 3942, 9,
                                   nmax);
    }

    for (size_t m = 0; m < 162; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
