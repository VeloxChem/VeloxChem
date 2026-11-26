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


#include "SimdThreeCenterElectronRepulsionGeom100RecPDF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
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
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_pdf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_pdf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 4492, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 315 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 4492, 2257, 1052, dimensions);

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
                                                        5, 6, 7}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 14, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 17, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 20, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 23, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 7, 8,
                                                                       14, 17, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 8, 9,
                                                                       17, 20, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 9, 10,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 10, 11,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 11, 12,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 14, 17,
                                                                       32, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 17, 20,
                                                                       38, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 20, 23,
                                                                       44, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 23, 26,
                                                                       50, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 32, 38,
                                                                       62, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 117, 0, 3, 38, 44,
                                                                       72, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 44, 50,
                                                                       82, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 147, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 150, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 153, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 156, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 159, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 162, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 165, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 168, 3, 7, 14,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 177, 3, 8, 17,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 186, 3, 9, 20,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 195, 3, 10, 23,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 204, 3, 11, 26,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 213, 3, 12, 29,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 222, 3, 14, 32,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 240, 3, 17, 38,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 258, 3, 20, 44,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 276, 3, 23, 50,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 294, 3, 26, 56,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 312, 3, 32, 62,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 342, 3, 38, 72,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 372, 3, 44, 82,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 402, 3, 50, 92,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 432, 3, 62, 102,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 477, 3, 72, 117,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 522, 3, 82, 132,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 567, 3, 7, 8, 153,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 573, 3, 8, 9, 156,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 579, 3, 9, 10,
                                                                       159, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 585, 3, 10, 11,
                                                                       162, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 591, 3, 11, 12,
                                                                       165, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 597, 0, 3, 567,
                                                                       153, 573, 186, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 615, 0, 3, 573,
                                                                       156, 579, 195, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 633, 0, 3, 579,
                                                                       159, 585, 204, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 651, 0, 3, 585,
                                                                       162, 591, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 669, 0, 3, 597,
                                                                       186, 615, 32, 38, 258,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 705, 0, 3, 615,
                                                                       195, 633, 38, 44, 276,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 741, 0, 3, 633,
                                                                       204, 651, 44, 50, 294,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 777, 0, 3, 669,
                                                                       258, 705, 62, 72, 372,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 837, 0, 3, 705,
                                                                       276, 741, 72, 82, 402,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 897, 0, 3, 777,
                                                                       372, 837, 102, 117, 522,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 987, 3, 147, 150,
                                                                       567, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 997, 3, 150, 153,
                                                                       573, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1007, 3, 153, 156,
                                                                       579, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1017, 3, 156, 159,
                                                                       585, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1027, 3, 159, 162,
                                                                       591, ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1037, 0, 3, 987,
                                                                       567, 997, 168, 177, 597,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1067, 0, 3, 997,
                                                                       573, 1007, 177, 186, 615,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1097, 0, 3, 1007,
                                                                       579, 1017, 186, 195, 633,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1127, 0, 3, 1017,
                                                                       585, 1027, 195, 204, 651,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1157, 0, 3, 1037,
                                                                       597, 1067, 222, 240, 669,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1217, 0, 3, 1067,
                                                                       615, 1097, 240, 258, 705,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1277, 0, 3, 1097,
                                                                       633, 1127, 258, 276, 741,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 1337, 0, 3, 1157,
                                                                       669, 1217, 312, 342, 777,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 1437, 0, 3, 1217,
                                                                       705, 1277, 342, 372, 837,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 1537, 0, 3, 1337,
                                                                       777, 1437, 432, 477, 897,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_p_x(buffer, 1687, 987, 1157, 1, 10, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 1717, 987, 1157, 1, 10, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 1747, 987, 1157, 1, 10, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 1777, 1037, 1337, 1, 10, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 1837, 1037, 1337, 1, 10, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 1897, 1037, 1337, 1, 10, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 1957, 1157, 1537, 1, 10, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 2057, 1157, 1537, 1, 10, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 2157, 1157, 1537, 1, 10, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 2257, 1687, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2308, 1717, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2359, 1747, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2410, 1037, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2461, 1777, 60, ncols);

                    simdfunc::contract_primitives(buffer, 2563, 1837, 60, ncols);

                    simdfunc::contract_primitives(buffer, 2665, 1897, 60, ncols);

                    simdfunc::contract_primitives(buffer, 2767, 1157, 60, ncols);

                    simdfunc::contract_primitives(buffer, 2869, 1957, 100, ncols);

                    simdfunc::contract_primitives(buffer, 3039, 2057, 100, ncols);

                    simdfunc::contract_primitives(buffer, 3209, 2157, 100, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 2287, 2257, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2338, 2308, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2389, 2359, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2440, 2410, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2521, 2461, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2623, 2563, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2725, 2665, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2827, 2767, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2969, 2869, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 3139, 3039, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 3309, 3209, 10, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 3379, 2287, 2440,
                                                       2521, 7, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 3442, 2338, 2440,
                                                       2623, 7, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 3505, 2389, 2440,
                                                       2725, 7, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 3568, 2440, 2827, 7, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 3631, 2521, 2827,
                                                       2969, 7, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 3757, 2623, 2827,
                                                       3139, 7, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 3883, 2725, 2827,
                                                       3309, 7, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 4009, 3379, 3568,
                                                       3631, 7, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 4135, 3442, 3568,
                                                       3757, 7, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 4261, 3505, 3568,
                                                       3883, 7, nmax);

        simdtrf::transform_d_inner(buffer, 4387, 4009, 3, 7, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 4387, 35, nmax);

        simdtrf::transform_d_inner(buffer, 4387, 4135, 3, 7, nmax);

        simdtrf::transform_p_outer(values + 105 * nvalues + n * npairs, nvalues, buffer, 4387,
                                   35, nmax);

        simdtrf::transform_d_inner(buffer, 4387, 4261, 3, 7, nmax);

        simdtrf::transform_p_outer(values + 210 * nvalues + n * npairs, nvalues, buffer, 4387,
                                   35, nmax);
    }

    for (size_t m = 0; m < 315; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
