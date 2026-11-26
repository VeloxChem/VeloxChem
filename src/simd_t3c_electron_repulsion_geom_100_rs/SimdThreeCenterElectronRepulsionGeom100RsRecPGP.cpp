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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecPGP.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
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

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_pgp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_pgp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 11555, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 486 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 11555, 1922, 2325, dimensions);

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
                                                            4, 5, 6, 7}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 14, 3, {1, 2, 3, 4,
                                                        5, 6, 7}, ncols, fj, i * nprim_b + j,
                                                        fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 428, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 431, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 434, 3, 7, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 443, 3, 15, 40,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 452, 3, 22, 58,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 470, 3, 40, 88,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 488, 3, 58, 118,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 518, 3, 88, 158,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 548, 3, 118, 198,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 593, 3, 158, 243,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 638, 3, 198, 288,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 701, 3, 243, 330,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 764, 3, 288, 372,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 848, 3, 330, 400,
                                                                       ncols, p, q);

                    simdgeo::geom_p_x(buffer, 932, 428, 452, 1, 3, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 941, 428, 452, 1, 3, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 950, 428, 452, 1, 3, ncols, alpha);

                    simdgeo::geom_p_x(buffer, 959, 431, 470, 1, 3, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 968, 431, 470, 1, 3, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 977, 431, 470, 1, 3, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 986, 434, 488, 1, 3, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 1004, 434, 488, 1, 3, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 1022, 434, 488, 1, 3, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 1040, 443, 518, 1, 3, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 1058, 443, 518, 1, 3, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 1076, 443, 518, 1, 3, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 1094, 452, 548, 1, 3, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 1124, 452, 548, 1, 3, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 1154, 452, 548, 1, 3, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 1184, 470, 593, 1, 3, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 1214, 470, 593, 1, 3, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 1244, 470, 593, 1, 3, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 1274, 488, 638, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 1319, 488, 638, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 1364, 488, 638, 1, 3, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 1409, 518, 701, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 1454, 518, 701, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 1499, 518, 701, 1, 3, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 1544, 548, 764, 1, 3, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 1607, 548, 764, 1, 3, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 1670, 548, 764, 1, 3, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 1733, 593, 848, 1, 3, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 1796, 593, 848, 1, 3, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 1859, 593, 848, 1, 3, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 1922, 932, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1940, 941, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1958, 950, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1976, 434, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1994, 959, 9, ncols);

                    simdfunc::contract_primitives(buffer, 2012, 968, 9, ncols);

                    simdfunc::contract_primitives(buffer, 2030, 977, 9, ncols);

                    simdfunc::contract_primitives(buffer, 2048, 443, 9, ncols);

                    simdfunc::contract_primitives(buffer, 2066, 986, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2102, 1004, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2138, 1022, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2174, 452, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2210, 1040, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2246, 1058, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2282, 1076, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2318, 470, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2354, 1094, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2414, 1124, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2474, 1154, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2534, 488, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2594, 1184, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2654, 1214, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2714, 1244, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2774, 518, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2834, 1274, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2924, 1319, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3014, 1364, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3104, 548, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3194, 1409, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3284, 1454, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3374, 1499, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3464, 593, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3554, 1544, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3680, 1607, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3806, 1670, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3932, 1733, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4058, 1796, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4184, 1859, 63, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 1931, 1922, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1949, 1940, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1967, 1958, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1985, 1976, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2003, 1994, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2021, 2012, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2039, 2030, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2057, 2048, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2084, 2066, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2120, 2102, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2156, 2138, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2192, 2174, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2228, 2210, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2264, 2246, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2300, 2282, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2336, 2318, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2384, 2354, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2444, 2414, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2504, 2474, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2564, 2534, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2624, 2594, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2684, 2654, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2744, 2714, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2804, 2774, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2879, 2834, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2969, 2924, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3059, 3014, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3149, 3104, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3239, 3194, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3329, 3284, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3419, 3374, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3509, 3464, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3617, 3554, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3743, 3680, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3869, 3806, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3995, 3932, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4121, 4058, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4247, 4184, 21, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 4310, 1931, 1985,
                                                       2084, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 4337, 1949, 1985,
                                                       2120, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 4364, 1967, 1985,
                                                       2156, 3, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 4391, 1985, 2192, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 4418, 2003, 2057,
                                                       2228, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 4445, 2021, 2057,
                                                       2264, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 4472, 2039, 2057,
                                                       2300, 3, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 4499, 2057, 2336, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 4526, 2084, 2192,
                                                       2384, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 4580, 2120, 2192,
                                                       2444, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 4634, 2156, 2192,
                                                       2504, 3, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 4688, 2192, 2564, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 4742, 2228, 2336,
                                                       2624, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 4796, 2264, 2336,
                                                       2684, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 4850, 2300, 2336,
                                                       2744, 3, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 4904, 2336, 2804, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 4958, 2384, 2564,
                                                       2879, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 5048, 2444, 2564,
                                                       2969, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 5138, 2504, 2564,
                                                       3059, 3, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 5228, 2564, 3149, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 5318, 2624, 2804,
                                                       3239, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 5408, 2684, 2804,
                                                       3329, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 5498, 2744, 2804,
                                                       3419, 3, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 5588, 2804, 3509, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 5678, 2879, 3149,
                                                       3617, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 5813, 2969, 3149,
                                                       3743, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 5948, 3059, 3149,
                                                       3869, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 6083, 3239, 3509,
                                                       3995, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 6218, 3329, 3509,
                                                       4121, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 6353, 3419, 3509,
                                                       4247, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 6488, 4310, 4391,
                                                       4526, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 6542, 4337, 4391,
                                                       4580, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 6596, 4364, 4391,
                                                       4634, 3, nmax);

        simdtrf::compute_hrr_pd_out_of_first(buffer, coordinates, 6650, 4391, 4688, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 6704, 4418, 4499,
                                                       4742, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 6758, 4445, 4499,
                                                       4796, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 6812, 4472, 4499,
                                                       4850, 3, nmax);

        simdtrf::compute_hrr_pd_out_of_first(buffer, coordinates, 6866, 4499, 4904, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 6920, 4526, 4688,
                                                       4958, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 7028, 4580, 4688,
                                                       5048, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 7136, 4634, 4688,
                                                       5138, 3, nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 7244, 4688, 5228, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 7352, 4742, 4904,
                                                       5318, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 7460, 4796, 4904,
                                                       5408, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 7568, 4850, 4904,
                                                       5498, 3, nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 7676, 4904, 5588, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 7784, 4958, 5228,
                                                       5678, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 7964, 5048, 5228,
                                                       5813, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 8144, 5138, 5228,
                                                       5948, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 8324, 5318, 5588,
                                                       6083, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 8504, 5408, 5588,
                                                       6218, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 8684, 5498, 5588,
                                                       6353, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pf_out_of_first(buffer, coordinates, 8864, 6488, 6650,
                                                       6920, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pf_out_of_first(buffer, coordinates, 8954, 6542, 6650,
                                                       7028, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pf_out_of_first(buffer, coordinates, 9044, 6596, 6650,
                                                       7136, 3, nmax);

        simdtrf::compute_hrr_pf_out_of_first(buffer, coordinates, 9134, 6650, 7244, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pf_out_of_first(buffer, coordinates, 9224, 6704, 6866,
                                                       7352, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pf_out_of_first(buffer, coordinates, 9314, 6758, 6866,
                                                       7460, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pf_out_of_first(buffer, coordinates, 9404, 6812, 6866,
                                                       7568, 3, nmax);

        simdtrf::compute_hrr_pf_out_of_first(buffer, coordinates, 9494, 6866, 7676, 3, nmax);

        simdtrf::compute_hrr_geom_100x_df_out_of_first(buffer, coordinates, 9584, 6920, 7244,
                                                       7784, 3, nmax);

        simdtrf::compute_hrr_geom_100y_df_out_of_first(buffer, coordinates, 9764, 7028, 7244,
                                                       7964, 3, nmax);

        simdtrf::compute_hrr_geom_100z_df_out_of_first(buffer, coordinates, 9944, 7136, 7244,
                                                       8144, 3, nmax);

        simdtrf::compute_hrr_geom_100x_df_out_of_first(buffer, coordinates, 10124, 7352, 7676,
                                                       8324, 3, nmax);

        simdtrf::compute_hrr_geom_100y_df_out_of_first(buffer, coordinates, 10304, 7460, 7676,
                                                       8504, 3, nmax);

        simdtrf::compute_hrr_geom_100z_df_out_of_first(buffer, coordinates, 10484, 7568, 7676,
                                                       8684, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pg_out_of_first(buffer, coordinates, 10664, 8864, 9134,
                                                       9584, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pg_out_of_first(buffer, coordinates, 10799, 8954, 9134,
                                                       9764, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pg_out_of_first(buffer, coordinates, 10934, 9044, 9134,
                                                       9944, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pg_out_of_first(buffer, coordinates, 11069, 9224, 9494,
                                                       10124, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pg_out_of_first(buffer, coordinates, 11204, 9314, 9494,
                                                       10304, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pg_out_of_first(buffer, coordinates, 11339, 9404, 9494,
                                                       10484, 3, nmax);

        simdtrf::transform_g_inner(buffer, 11474, 11069, 3, 3, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 11474, 27, nmax);

        simdtrf::transform_g_inner(buffer, 11474, 11204, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 81 * nvalues + n * npairs, nvalues, buffer, 11474,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 11474, 11339, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 162 * nvalues + n * npairs, nvalues, buffer, 11474,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 11474, 10664, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 243 * nvalues + n * npairs, nvalues, buffer, 11474,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 11474, 10799, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 324 * nvalues + n * npairs, nvalues, buffer, 11474,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 11474, 10934, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 405 * nvalues + n * npairs, nvalues, buffer, 11474,
                                   27, nmax);
    }

    for (size_t m = 0; m < 486; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
