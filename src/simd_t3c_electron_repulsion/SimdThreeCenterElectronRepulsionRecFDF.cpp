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


#include "SimdThreeCenterElectronRepulsionRecFDF.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferFD.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_fdf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_fdf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto dimensions = simdfunc::make_column_dimensions(
        a_function, b_function, c_function, npairs, coordinates,
        screenfunc::three_center_electron_repulsion_primitive_bound,
        threshold / static_cast<double>(nprims));

    const auto nmax = simdfunc::prepare_buffer(buffer, 5053, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 245 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 5053, 2976, 635, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fa = -b_exps[j] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pa(buffer, coordinates, 0, nmax, fa);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 3, n, nmax, fc);

                simdfunc::compute_pair_exponent(buffer, coordinates, 6, nmax, mu);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 7, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8}, ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 16, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 19, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 22, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 8, 9,
                                                                       16, 19, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 9, 10,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 10, 11,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 11, 12,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 12, 13,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 13, 14,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 16, 19,
                                                                       37, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 19, 22,
                                                                       43, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 93, 0, 3, 22, 25,
                                                                       49, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 25, 28,
                                                                       55, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 28, 31,
                                                                       61, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 123, 0, 3, 37, 43,
                                                                       73, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 43, 49,
                                                                       83, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 153, 0, 3, 49, 55,
                                                                       93, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 168, 0, 3, 55, 61,
                                                                       103, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 73, 83,
                                                                       123, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 204, 0, 3, 83, 93,
                                                                       138, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 225, 0, 3, 93,
                                                                       103, 153, 168, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 246, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 249, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 252, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 255, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 258, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 261, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 264, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 267, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 270, 3, 8, 16,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 279, 3, 9, 19,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 288, 3, 10, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 297, 3, 11, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 306, 3, 12, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 315, 3, 13, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 324, 3, 14, 34,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 333, 3, 16, 37,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 351, 3, 19, 43,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 369, 3, 22, 49,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 387, 3, 25, 55,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 405, 3, 28, 61,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 423, 3, 31, 67,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 441, 3, 37, 73,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 471, 3, 43, 83,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 501, 3, 49, 93,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 531, 3, 55, 103,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 561, 3, 61, 113,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 591, 3, 73, 123,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 636, 3, 83, 138,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 681, 3, 93, 153,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 726, 3, 103, 168,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 771, 3, 123, 183,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 834, 3, 138, 204,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 897, 3, 153, 225,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 960, 3, 8, 9, 252,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 966, 3, 9, 10,
                                                                       255, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 972, 3, 10, 11,
                                                                       258, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 978, 3, 11, 12,
                                                                       261, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 984, 3, 12, 13,
                                                                       264, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 990, 3, 13, 14,
                                                                       267, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 996, 0, 3, 960,
                                                                       252, 966, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1014, 0, 3, 966,
                                                                       255, 972, 297, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1032, 0, 3, 972,
                                                                       258, 978, 306, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1050, 0, 3, 978,
                                                                       261, 984, 315, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1068, 0, 3, 984,
                                                                       264, 990, 324, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1086, 0, 3, 996,
                                                                       288, 1014, 37, 43, 369,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1122, 0, 3, 1014,
                                                                       297, 1032, 43, 49, 387,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1158, 0, 3, 1032,
                                                                       306, 1050, 49, 55, 405,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1194, 0, 3, 1050,
                                                                       315, 1068, 55, 61, 423,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1230, 0, 3, 1086,
                                                                       369, 1122, 73, 83, 501,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1290, 0, 3, 1122,
                                                                       387, 1158, 83, 93, 531,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1350, 0, 3, 1158,
                                                                       405, 1194, 93, 103, 561,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 1410, 0, 3, 1230,
                                                                       501, 1290, 123, 138, 681,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 1500, 0, 3, 1290,
                                                                       531, 1350, 138, 153, 726,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 1590, 0, 3, 1410,
                                                                       681, 1500, 183, 204, 897,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1716, 3, 246, 249,
                                                                       960, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1726, 3, 249, 252,
                                                                       966, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1736, 3, 252, 255,
                                                                       972, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1746, 3, 255, 258,
                                                                       978, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1756, 3, 258, 261,
                                                                       984, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1766, 3, 261, 264,
                                                                       990, ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1776, 0, 3, 1716,
                                                                       960, 1726, 270, 279, 996,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1806, 0, 3, 1726,
                                                                       966, 1736, 279, 288, 1014,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1836, 0, 3, 1736,
                                                                       972, 1746, 288, 297, 1032,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1866, 0, 3, 1746,
                                                                       978, 1756, 297, 306, 1050,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1896, 0, 3, 1756,
                                                                       984, 1766, 306, 315, 1068,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1926, 0, 3, 1776,
                                                                       996, 1806, 333, 351, 1086,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1986, 0, 3, 1806,
                                                                       1014, 1836, 351, 369,
                                                                       1122, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2046, 0, 3, 1836,
                                                                       1032, 1866, 369, 387,
                                                                       1158, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2106, 0, 3, 1866,
                                                                       1050, 1896, 387, 405,
                                                                       1194, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 2166, 0, 3, 1926,
                                                                       1086, 1986, 441, 471,
                                                                       1230, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 2266, 0, 3, 1986,
                                                                       1122, 2046, 471, 501,
                                                                       1290, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 2366, 0, 3, 2046,
                                                                       1158, 2106, 501, 531,
                                                                       1350, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 2466, 0, 3, 2166,
                                                                       1230, 2266, 591, 636,
                                                                       1410, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 2616, 0, 3, 2266,
                                                                       1290, 2366, 636, 681,
                                                                       1500, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 2766, 0, 3, 2466,
                                                                       1410, 2616, 771, 834,
                                                                       1590, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 2976, 2166, 100, ncols);

                    simdfunc::contract_primitives(buffer, 3146, 2466, 150, ncols);

                    simdfunc::contract_primitives(buffer, 3401, 2766, 210, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 3076, 2976, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 3296, 3146, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 3611, 3401, 21, 1, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 3758, 3076, 3296, 7, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 3968, 3296, 3611, 7, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 4283, 3758, 3968, 7, nmax);

        simdtrf::transform_d_inner(buffer, 4703, 4283, 10, 7, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 4703, 35, nmax);
    }

    for (size_t m = 0; m < 245; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
