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


#include "SimdElectronRepulsionGeom10RsRecFF.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdElectronRepulsionVrrRecDD.hpp"
#include "SimdElectronRepulsionVrrRecDF.hpp"
#include "SimdElectronRepulsionVrrRecDP.hpp"
#include "SimdElectronRepulsionVrrRecDS.hpp"
#include "SimdElectronRepulsionVrrRecFD.hpp"
#include "SimdElectronRepulsionVrrRecFF.hpp"
#include "SimdElectronRepulsionVrrRecFP.hpp"
#include "SimdElectronRepulsionVrrRecFS.hpp"
#include "SimdElectronRepulsionVrrRecGD.hpp"
#include "SimdElectronRepulsionVrrRecGF.hpp"
#include "SimdElectronRepulsionVrrRecGP.hpp"
#include "SimdElectronRepulsionVrrRecGS.hpp"
#include "SimdElectronRepulsionVrrRecPD.hpp"
#include "SimdElectronRepulsionVrrRecPF.hpp"
#include "SimdElectronRepulsionVrrRecPP.hpp"
#include "SimdElectronRepulsionVrrRecPS.hpp"
#include "SimdElectronRepulsionVrrRecSD.hpp"
#include "SimdElectronRepulsionVrrRecSF.hpp"
#include "SimdElectronRepulsionVrrRecSP.hpp"
#include "SimdGeometryF1.hpp"
#include "SimdTransformF.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

auto
compute_rs_geom_10_ff_electron_repulsion(double               *values,
                                         const size_t          nvalues,
                                         const CBasisFunction &bra,
                                         const CBasisFunction &ket,
                                         const CSimdMatrix    &coordinates,
                                         CSimdMatrix          &buffer,
                                         const double          omega) -> void
{
    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_geom_10_ff_electron_repulsion: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nmax = simdfunc::prepare_buffer(buffer, 4236, 3566, 600, nvalues);

    for (size_t i = 0; i < nprim_a; i++)
    {
        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto ncols = nvalues;

            const auto p = a_exps[i] + b_exps[j];

            const auto mu = a_exps[i] * b_exps[j] / p;

            const auto pi = mathconst::pi_value();

            const auto fj = 2.0 * a_norms[i] * b_norms[j] * pi * pi * std::sqrt(pi)
                            / (a_exps[i] * b_exps[j] * std::sqrt(p));

            const auto alpha = a_exps[i];

            const auto beta = b_exps[j];

            const auto fb = a_exps[i] / p;

            const auto fa = -b_exps[j] / p;

            simdfunc::compute_pa(buffer, coordinates, 0, ncols, fa);

            simdfunc::compute_pb(buffer, coordinates, 3, ncols, fb);

            simdfunc::compute_erf_boys_function(buffer, coordinates, 6, {1, 2, 3, 4, 5, 6, 7},
                                                ncols, fj, mu, omega);

            simdfunc::compute_boys_function(buffer, coordinates, 14, {1, 2, 3, 4, 5, 6, 7},
                                            ncols, fj, mu);

            compute_prim_ps_electron_repulsion_0(buffer, 22, 0, 8, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 25, 0, 9, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 28, 0, 10, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 31, 0, 11, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 34, 0, 12, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 37, 0, 13, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 40, 0, 16, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 43, 0, 17, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 46, 0, 18, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 49, 0, 19, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 52, 0, 20, ncols);

            compute_prim_ps_electron_repulsion_0(buffer, 55, 0, 21, ncols);

            compute_prim_ds_electron_repulsion_0(buffer, 58, 0, 7, 8, 25, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 64, 0, 8, 9, 28, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 70, 0, 9, 10, 31, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 76, 0, 10, 11, 34, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 82, 0, 11, 12, 37, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 88, 0, 15, 16, 43, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 94, 0, 16, 17, 46, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 100, 0, 17, 18, 49, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 106, 0, 18, 19, 52, ncols, alpha, beta,
                                                 p);

            compute_prim_ds_electron_repulsion_0(buffer, 112, 0, 19, 20, 55, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 118, 0, 22, 25, 64, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 128, 0, 25, 28, 70, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 138, 0, 28, 31, 76, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 148, 0, 31, 34, 82, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 158, 0, 40, 43, 94, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 168, 0, 43, 46, 100, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 178, 0, 46, 49, 106, ncols, alpha, beta,
                                                 p);

            compute_prim_fs_electron_repulsion_0(buffer, 188, 0, 49, 52, 112, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 198, 0, 58, 64, 128, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 213, 0, 64, 70, 138, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 228, 0, 70, 76, 148, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 243, 0, 88, 94, 168, ncols, alpha, beta,
                                                 p);

            compute_prim_gs_electron_repulsion_0(buffer, 258, 0, 94, 100, 178, ncols, alpha,
                                                 beta, p);

            compute_prim_gs_electron_repulsion_0(buffer, 273, 0, 100, 106, 188, ncols, alpha,
                                                 beta, p);

            compute_prim_sp_electron_repulsion_0(buffer, 288, 3, 9, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 291, 3, 10, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 294, 3, 11, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 297, 3, 12, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 300, 3, 13, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 303, 3, 17, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 306, 3, 18, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 309, 3, 19, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 312, 3, 20, ncols);

            compute_prim_sp_electron_repulsion_0(buffer, 315, 3, 21, ncols);

            compute_prim_pp_electron_repulsion_0(buffer, 318, 3, 8, 25, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 327, 3, 9, 28, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 336, 3, 10, 31, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 345, 3, 11, 34, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 354, 3, 12, 37, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 363, 3, 16, 43, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 372, 3, 17, 46, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 381, 3, 18, 49, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 390, 3, 19, 52, ncols, p);

            compute_prim_pp_electron_repulsion_0(buffer, 399, 3, 20, 55, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 408, 0, 3, 22, 318, 58, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 426, 0, 3, 25, 327, 64, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 444, 0, 3, 28, 336, 70, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 462, 0, 3, 31, 345, 76, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 480, 0, 3, 34, 354, 82, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 498, 0, 3, 40, 363, 88, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 516, 0, 3, 43, 372, 94, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 534, 0, 3, 46, 381, 100, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 552, 0, 3, 49, 390, 106, ncols, p);

            compute_prim_dp_electron_repulsion_0(buffer, 570, 0, 3, 52, 399, 112, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 588, 0, 3, 64, 444, 128, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 618, 0, 3, 70, 462, 138, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 648, 0, 3, 76, 480, 148, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 678, 0, 3, 94, 534, 168, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 708, 0, 3, 100, 552, 178, ncols, p);

            compute_prim_fp_electron_repulsion_0(buffer, 738, 0, 3, 106, 570, 188, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 768, 0, 3, 118, 588, 198, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 813, 0, 3, 128, 618, 213, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 858, 0, 3, 138, 648, 228, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 903, 0, 3, 158, 678, 243, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 948, 0, 3, 168, 708, 258, ncols, p);

            compute_prim_gp_electron_repulsion_0(buffer, 993, 0, 3, 178, 738, 273, ncols, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1038, 3, 8, 9, 291, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1044, 3, 9, 10, 294, ncols, alpha, beta,
                                                 p);

            compute_prim_sd_electron_repulsion_0(buffer, 1050, 3, 10, 11, 297, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1056, 3, 11, 12, 300, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1062, 3, 16, 17, 306, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1068, 3, 17, 18, 309, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1074, 3, 18, 19, 312, ncols, alpha,
                                                 beta, p);

            compute_prim_sd_electron_repulsion_0(buffer, 1080, 3, 19, 20, 315, ncols, alpha,
                                                 beta, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1086, 0, 3, 288, 1038, 327, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1104, 0, 3, 291, 1044, 336, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1122, 0, 3, 294, 1050, 345, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1140, 0, 3, 297, 1056, 354, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1158, 0, 3, 303, 1062, 372, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1176, 0, 3, 306, 1068, 381, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1194, 0, 3, 309, 1074, 390, ncols, p);

            compute_prim_pd_electron_repulsion_0(buffer, 1212, 0, 3, 312, 1080, 399, ncols, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1230, 0, 3, 327, 1104, 58, 64, 444,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1266, 0, 3, 336, 1122, 64, 70, 462,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1302, 0, 3, 345, 1140, 70, 76, 480,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1338, 0, 3, 372, 1176, 88, 94, 534,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1374, 0, 3, 381, 1194, 94, 100, 552,
                                                 ncols, alpha, beta, p);

            compute_prim_dd_electron_repulsion_0(buffer, 1410, 0, 3, 390, 1212, 100, 106, 570,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1446, 0, 3, 444, 1266, 118, 128, 618,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1506, 0, 3, 462, 1302, 128, 138, 648,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1566, 0, 3, 534, 1374, 158, 168, 708,
                                                 ncols, alpha, beta, p);

            compute_prim_fd_electron_repulsion_0(buffer, 1626, 0, 3, 552, 1410, 168, 178, 738,
                                                 ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1686, 0, 3, 1230, 1266, 618, 1506, 198,
                                                 213, 858, ncols, alpha, beta, p);

            compute_prim_gd_electron_repulsion_0(buffer, 1776, 0, 3, 1338, 1374, 708, 1626, 243,
                                                 258, 993, ncols, alpha, beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1866, 3, 288, 291, 1044, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1876, 3, 291, 294, 1050, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1886, 3, 294, 297, 1056, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1896, 3, 303, 306, 1068, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1906, 3, 306, 309, 1074, ncols, alpha,
                                                 beta, p);

            compute_prim_sf_electron_repulsion_0(buffer, 1916, 3, 309, 312, 1080, ncols, alpha,
                                                 beta, p);

            compute_prim_pf_electron_repulsion_0(buffer, 1926, 0, 3, 1038, 1866, 1104, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 1956, 0, 3, 1044, 1876, 1122, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 1986, 0, 3, 1050, 1886, 1140, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2016, 0, 3, 1062, 1896, 1176, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2046, 0, 3, 1068, 1906, 1194, ncols,
                                                 p);

            compute_prim_pf_electron_repulsion_0(buffer, 2076, 0, 3, 1074, 1916, 1212, ncols,
                                                 p);

            compute_prim_df_electron_repulsion_0(buffer, 2106, 0, 3, 1086, 1926, 408, 426, 1230,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2166, 0, 3, 1104, 1956, 426, 444, 1266,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2226, 0, 3, 1122, 1986, 444, 462, 1302,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2286, 0, 3, 1158, 2016, 498, 516, 1338,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2346, 0, 3, 1176, 2046, 516, 534, 1374,
                                                 ncols, alpha, beta, p);

            compute_prim_df_electron_repulsion_0(buffer, 2406, 0, 3, 1194, 2076, 534, 552, 1410,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2466, 0, 3, 1266, 2226, 588, 618, 1506,
                                                 ncols, alpha, beta, p);

            compute_prim_ff_electron_repulsion_0(buffer, 2566, 0, 3, 1374, 2406, 678, 708, 1626,
                                                 ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 2666, 0, 3, 2106, 2166, 1446, 2466, 768,
                                                 813, 1686, ncols, alpha, beta, p);

            compute_prim_gf_electron_repulsion_0(buffer, 2816, 0, 3, 2286, 2346, 1566, 2566, 903,
                                                 948, 1776, ncols, alpha, beta, p);

            simdgeo::geom_f_x(buffer, 2966, 2286, 2816, 1, 10, ncols, alpha);

            simdgeo::geom_f_y(buffer, 3066, 2286, 2816, 1, 10, ncols, alpha);

            simdgeo::geom_f_z(buffer, 3166, 2286, 2816, 1, 10, ncols, alpha);

            simdgeo::geom_f_x(buffer, 3266, 2106, 2666, 1, 10, ncols, alpha);

            simdgeo::geom_f_y(buffer, 3366, 2106, 2666, 1, 10, ncols, alpha);

            simdgeo::geom_f_z(buffer, 3466, 2106, 2666, 1, 10, ncols, alpha);

            simdfunc::contract_primitives(buffer, 3566, 3266, 300, ncols);

            simdfunc::contract_primitives(buffer, 3866, 2966, 300, ncols);
        }
    }

    simdtrf::transform_f_inner(buffer, 4166, 3866, 10, 1, nmax);

    simdtrf::transform_f_outer(values, nvalues, buffer, 4166, 7, nmax);

    simdtrf::transform_f_inner(buffer, 4166, 3966, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 49 * nvalues, nvalues, buffer, 4166, 7, nmax);

    simdtrf::transform_f_inner(buffer, 4166, 4066, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 98 * nvalues, nvalues, buffer, 4166, 7, nmax);

    simdtrf::transform_f_inner(buffer, 4166, 3566, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 147 * nvalues, nvalues, buffer, 4166, 7, nmax);

    simdtrf::transform_f_inner(buffer, 4166, 3666, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 196 * nvalues, nvalues, buffer, 4166, 7, nmax);

    simdtrf::transform_f_inner(buffer, 4166, 3766, 10, 1, nmax);

    simdtrf::transform_f_outer(values + 245 * nvalues, nvalues, buffer, 4166, 7, nmax);
}

}  // namespace simdt2ceri
