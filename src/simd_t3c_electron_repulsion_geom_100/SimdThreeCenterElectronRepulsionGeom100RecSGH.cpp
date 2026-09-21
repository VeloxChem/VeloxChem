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


#include "SimdThreeCenterElectronRepulsionGeom100RecSGH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGS.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_sgh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_sgh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 26407, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 297 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 26407, 25297, 945, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto alpha = a_exps[i];

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
                                                        5, 6, 7, 8, 9, 10}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 20, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 47, 3, 6, 10, 11,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 53, 3, 6, 11, 12,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 59, 3, 6, 12, 13,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 65, 3, 6, 13, 14,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 71, 3, 6, 14, 15,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 77, 3, 6, 15, 16,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 83, 3, 6, 16, 17,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 89, 3, 6, 17, 18,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 95, 3, 6, 20, 23,
                                                                       47, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 105, 3, 6, 23, 26,
                                                                       53, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 115, 3, 6, 26, 29,
                                                                       59, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 125, 3, 6, 29, 32,
                                                                       65, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 135, 3, 6, 32, 35,
                                                                       71, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 145, 3, 6, 35, 38,
                                                                       77, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 155, 3, 6, 38, 41,
                                                                       83, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 165, 3, 6, 47, 53,
                                                                       95, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 180, 3, 6, 53, 59,
                                                                       105, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 195, 3, 6, 59, 65,
                                                                       115, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 210, 3, 6, 65, 71,
                                                                       125, 135, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 225, 3, 6, 71, 77,
                                                                       135, 145, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 240, 3, 6, 77, 83,
                                                                       145, 155, ncols, gamma, p,
                                                                       q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 255, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 258, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 261, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 264, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 267, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 270, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 273, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 276, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 279, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 282, 0, 6, 10, 11,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 291, 0, 6, 11, 12,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 300, 0, 6, 12, 13,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 309, 0, 6, 13, 14,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 318, 0, 6, 14, 15,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 327, 0, 6, 15, 16,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 336, 0, 6, 16, 17,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 345, 0, 6, 17, 18,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 354, 0, 3, 6, 20,
                                                                       23, 47, 53, 282, 291,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 372, 0, 3, 6, 23,
                                                                       26, 53, 59, 291, 300,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 390, 0, 3, 6, 26,
                                                                       29, 59, 65, 300, 309,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 408, 0, 3, 6, 29,
                                                                       32, 65, 71, 309, 318,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 426, 0, 3, 6, 32,
                                                                       35, 71, 77, 318, 327,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 444, 0, 3, 6, 35,
                                                                       38, 77, 83, 327, 336,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 462, 0, 3, 6, 38,
                                                                       41, 83, 89, 336, 345,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 480, 0, 3, 6, 47,
                                                                       53, 95, 105, 354, 372,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 510, 0, 3, 6, 53,
                                                                       59, 105, 115, 372, 390,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 540, 0, 3, 6, 59,
                                                                       65, 115, 125, 390, 408,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 570, 0, 3, 6, 65,
                                                                       71, 125, 135, 408, 426,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 600, 0, 3, 6, 71,
                                                                       77, 135, 145, 426, 444,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 630, 0, 3, 6, 77,
                                                                       83, 145, 155, 444, 462,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 660, 0, 3, 6, 95,
                                                                       105, 165, 180, 480, 510,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 705, 0, 3, 6, 105,
                                                                       115, 180, 195, 510, 540,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 750, 0, 3, 6, 115,
                                                                       125, 195, 210, 540, 570,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 795, 0, 3, 6, 125,
                                                                       135, 210, 225, 570, 600,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 840, 0, 3, 6, 135,
                                                                       145, 225, 240, 600, 630,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 885, 6, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 888, 6, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 891, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 894, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 897, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 900, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 903, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 906, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 909, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 912, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 915, 6, 12, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 924, 6, 13, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 933, 6, 14, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 942, 6, 15, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 951, 6, 16, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 960, 6, 17, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 969, 6, 18, 44,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 978, 6, 20, 47,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 996, 6, 23, 53,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1014, 6, 26, 59,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1032, 6, 29, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1050, 6, 32, 71,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1068, 6, 35, 77,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1086, 6, 38, 83,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1104, 6, 41, 89,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1122, 6, 47, 95,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1152, 6, 53, 105,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1182, 6, 59, 115,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1212, 6, 65, 125,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1242, 6, 71, 135,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1272, 6, 77, 145,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1302, 6, 83, 155,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1332, 6, 95, 165,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1377, 6, 105, 180,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1422, 6, 115, 195,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1467, 6, 125, 210,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1512, 6, 135, 225,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1557, 6, 145, 240,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1602, 6, 10, 255,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1611, 6, 11, 258,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1620, 6, 12, 261,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1629, 6, 13, 264,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1638, 6, 14, 267,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1647, 6, 15, 270,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1656, 6, 16, 273,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1665, 6, 17, 276,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1674, 6, 18, 279,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1683, 6, 20, 255,
                                                                       282, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1710, 6, 23, 258,
                                                                       291, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1737, 6, 26, 261,
                                                                       300, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1764, 6, 29, 264,
                                                                       309, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1791, 6, 32, 267,
                                                                       318, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1818, 6, 35, 270,
                                                                       327, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1845, 6, 38, 273,
                                                                       336, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1872, 6, 41, 276,
                                                                       345, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1899, 3, 6, 47,
                                                                       1683, 282, 1710, 354,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1953, 3, 6, 53,
                                                                       1710, 291, 1737, 372,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2007, 3, 6, 59,
                                                                       1737, 300, 1764, 390,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2061, 3, 6, 65,
                                                                       1764, 309, 1791, 408,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2115, 3, 6, 71,
                                                                       1791, 318, 1818, 426,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2169, 3, 6, 77,
                                                                       1818, 327, 1845, 444,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2223, 3, 6, 83,
                                                                       1845, 336, 1872, 462,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2277, 3, 6, 95,
                                                                       1899, 354, 1953, 480,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2367, 3, 6, 105,
                                                                       1953, 372, 2007, 510,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2457, 3, 6, 115,
                                                                       2007, 390, 2061, 540,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2547, 3, 6, 125,
                                                                       2061, 408, 2115, 570,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2637, 3, 6, 135,
                                                                       2115, 426, 2169, 600,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2727, 3, 6, 145,
                                                                       2169, 444, 2223, 630,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 2817, 3, 6, 165,
                                                                       2277, 480, 2367, 660,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 2952, 3, 6, 180,
                                                                       2367, 510, 2457, 705,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 3087, 3, 6, 195,
                                                                       2457, 540, 2547, 750,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 3222, 3, 6, 210,
                                                                       2547, 570, 2637, 795,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 3357, 3, 6, 225,
                                                                       2637, 600, 2727, 840,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3492, 6, 10, 11,
                                                                       891, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3498, 6, 11, 12,
                                                                       894, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3504, 6, 12, 13,
                                                                       897, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3510, 6, 13, 14,
                                                                       900, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3516, 6, 14, 15,
                                                                       903, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3522, 6, 15, 16,
                                                                       906, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3528, 6, 16, 17,
                                                                       909, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3534, 6, 17, 18,
                                                                       912, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3540, 3, 6, 3492,
                                                                       891, 3498, 915, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3558, 3, 6, 3498,
                                                                       894, 3504, 924, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3576, 3, 6, 3504,
                                                                       897, 3510, 933, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3594, 3, 6, 3510,
                                                                       900, 3516, 942, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3612, 3, 6, 3516,
                                                                       903, 3522, 951, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3630, 3, 6, 3522,
                                                                       906, 3528, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3648, 3, 6, 3528,
                                                                       909, 3534, 969, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3666, 3, 6, 3540,
                                                                       915, 3558, 47, 53, 1014,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3702, 3, 6, 3558,
                                                                       924, 3576, 53, 59, 1032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3738, 3, 6, 3576,
                                                                       933, 3594, 59, 65, 1050,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3774, 3, 6, 3594,
                                                                       942, 3612, 65, 71, 1068,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3810, 3, 6, 3612,
                                                                       951, 3630, 71, 77, 1086,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3846, 3, 6, 3630,
                                                                       960, 3648, 77, 83, 1104,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3882, 3, 6, 3666,
                                                                       1014, 3702, 95, 105, 1182,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3942, 3, 6, 3702,
                                                                       1032, 3738, 105, 115,
                                                                       1212, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4002, 3, 6, 3738,
                                                                       1050, 3774, 115, 125,
                                                                       1242, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4062, 3, 6, 3774,
                                                                       1068, 3810, 125, 135,
                                                                       1272, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4122, 3, 6, 3810,
                                                                       1086, 3846, 135, 145,
                                                                       1302, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4182, 3, 6, 3882,
                                                                       1182, 3942, 165, 180,
                                                                       1422, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4272, 3, 6, 3942,
                                                                       1212, 4002, 180, 195,
                                                                       1467, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4362, 3, 6, 4002,
                                                                       1242, 4062, 195, 210,
                                                                       1512, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4452, 3, 6, 4062,
                                                                       1272, 4122, 210, 225,
                                                                       1557, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4542, 0, 6, 3492,
                                                                       891, 3498, 1620, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4560, 0, 6, 3498,
                                                                       894, 3504, 1629, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4578, 0, 6, 3504,
                                                                       897, 3510, 1638, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4596, 0, 6, 3510,
                                                                       900, 3516, 1647, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4614, 0, 6, 3516,
                                                                       903, 3522, 1656, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4632, 0, 6, 3522,
                                                                       906, 3528, 1665, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4650, 0, 6, 3528,
                                                                       909, 3534, 1674, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4668, 0, 3, 6,
                                                                       3540, 915, 3558, 4542,
                                                                       1620, 4560, 282, 291,
                                                                       1737, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4722, 0, 3, 6,
                                                                       3558, 924, 3576, 4560,
                                                                       1629, 4578, 291, 300,
                                                                       1764, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4776, 0, 3, 6,
                                                                       3576, 933, 3594, 4578,
                                                                       1638, 4596, 300, 309,
                                                                       1791, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4830, 0, 3, 6,
                                                                       3594, 942, 3612, 4596,
                                                                       1647, 4614, 309, 318,
                                                                       1818, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4884, 0, 3, 6,
                                                                       3612, 951, 3630, 4614,
                                                                       1656, 4632, 318, 327,
                                                                       1845, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4938, 0, 3, 6,
                                                                       3630, 960, 3648, 4632,
                                                                       1665, 4650, 327, 336,
                                                                       1872, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4992, 0, 3, 6,
                                                                       3666, 1014, 3702, 4668,
                                                                       1737, 4722, 354, 372,
                                                                       2007, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5100, 0, 3, 6,
                                                                       3702, 1032, 3738, 4722,
                                                                       1764, 4776, 372, 390,
                                                                       2061, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5208, 0, 3, 6,
                                                                       3738, 1050, 3774, 4776,
                                                                       1791, 4830, 390, 408,
                                                                       2115, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5316, 0, 3, 6,
                                                                       3774, 1068, 3810, 4830,
                                                                       1818, 4884, 408, 426,
                                                                       2169, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5424, 0, 3, 6,
                                                                       3810, 1086, 3846, 4884,
                                                                       1845, 4938, 426, 444,
                                                                       2223, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 5532, 0, 3, 6,
                                                                       3882, 1182, 3942, 4668,
                                                                       4722, 4992, 2007, 5100,
                                                                       480, 510, 2457, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 5712, 0, 3, 6,
                                                                       3942, 1212, 4002, 4722,
                                                                       4776, 5100, 2061, 5208,
                                                                       510, 540, 2547, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 5892, 0, 3, 6,
                                                                       4002, 1242, 4062, 4776,
                                                                       4830, 5208, 2115, 5316,
                                                                       540, 570, 2637, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 6072, 0, 3, 6,
                                                                       4062, 1272, 4122, 4830,
                                                                       4884, 5316, 2169, 5424,
                                                                       570, 600, 2727, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 6252, 0, 3, 6,
                                                                       4182, 1422, 4272, 4992,
                                                                       5100, 5532, 2457, 5712,
                                                                       660, 705, 3087, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 6522, 0, 3, 6,
                                                                       4272, 1467, 4362, 5100,
                                                                       5208, 5712, 2547, 5892,
                                                                       705, 750, 3222, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 6792, 0, 3, 6,
                                                                       4362, 1512, 4452, 5208,
                                                                       5316, 5892, 2637, 6072,
                                                                       750, 795, 3357, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7062, 6, 885, 888,
                                                                       3492, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7072, 6, 888, 891,
                                                                       3498, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7082, 6, 891, 894,
                                                                       3504, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7092, 6, 894, 897,
                                                                       3510, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7102, 6, 897, 900,
                                                                       3516, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7112, 6, 900, 903,
                                                                       3522, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7122, 6, 903, 906,
                                                                       3528, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7132, 6, 906, 909,
                                                                       3534, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7142, 3, 6, 7062,
                                                                       3492, 7072, 3540, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7172, 3, 6, 7072,
                                                                       3498, 7082, 3558, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7202, 3, 6, 7082,
                                                                       3504, 7092, 3576, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7232, 3, 6, 7092,
                                                                       3510, 7102, 3594, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7262, 3, 6, 7102,
                                                                       3516, 7112, 3612, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7292, 3, 6, 7112,
                                                                       3522, 7122, 3630, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7322, 3, 6, 7122,
                                                                       3528, 7132, 3648, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7352, 3, 6, 7142,
                                                                       3540, 7172, 978, 996,
                                                                       3666, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7412, 3, 6, 7172,
                                                                       3558, 7202, 996, 1014,
                                                                       3702, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7472, 3, 6, 7202,
                                                                       3576, 7232, 1014, 1032,
                                                                       3738, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7532, 3, 6, 7232,
                                                                       3594, 7262, 1032, 1050,
                                                                       3774, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7592, 3, 6, 7262,
                                                                       3612, 7292, 1050, 1068,
                                                                       3810, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7652, 3, 6, 7292,
                                                                       3630, 7322, 1068, 1086,
                                                                       3846, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7712, 3, 6, 7352,
                                                                       3666, 7412, 1122, 1152,
                                                                       3882, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7812, 3, 6, 7412,
                                                                       3702, 7472, 1152, 1182,
                                                                       3942, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7912, 3, 6, 7472,
                                                                       3738, 7532, 1182, 1212,
                                                                       4002, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 8012, 3, 6, 7532,
                                                                       3774, 7592, 1212, 1242,
                                                                       4062, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 8112, 3, 6, 7592,
                                                                       3810, 7652, 1242, 1272,
                                                                       4122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8212, 3, 6, 7712,
                                                                       3882, 7812, 1332, 1377,
                                                                       4182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8362, 3, 6, 7812,
                                                                       3942, 7912, 1377, 1422,
                                                                       4272, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8512, 3, 6, 7912,
                                                                       4002, 8012, 1422, 1467,
                                                                       4362, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 8662, 3, 6, 8012,
                                                                       4062, 8112, 1467, 1512,
                                                                       4452, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8812, 0, 6, 7062,
                                                                       3492, 7072, 1602, 1611,
                                                                       4542, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8842, 0, 6, 7072,
                                                                       3498, 7082, 1611, 1620,
                                                                       4560, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8872, 0, 6, 7082,
                                                                       3504, 7092, 1620, 1629,
                                                                       4578, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8902, 0, 6, 7092,
                                                                       3510, 7102, 1629, 1638,
                                                                       4596, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8932, 0, 6, 7102,
                                                                       3516, 7112, 1638, 1647,
                                                                       4614, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8962, 0, 6, 7112,
                                                                       3522, 7122, 1647, 1656,
                                                                       4632, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8992, 0, 6, 7122,
                                                                       3528, 7132, 1656, 1665,
                                                                       4650, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9022, 0, 3, 6,
                                                                       7142, 3540, 7172, 8812,
                                                                       4542, 8842, 1683, 1710,
                                                                       4668, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9112, 0, 3, 6,
                                                                       7172, 3558, 7202, 8842,
                                                                       4560, 8872, 1710, 1737,
                                                                       4722, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9202, 0, 3, 6,
                                                                       7202, 3576, 7232, 8872,
                                                                       4578, 8902, 1737, 1764,
                                                                       4776, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9292, 0, 3, 6,
                                                                       7232, 3594, 7262, 8902,
                                                                       4596, 8932, 1764, 1791,
                                                                       4830, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9382, 0, 3, 6,
                                                                       7262, 3612, 7292, 8932,
                                                                       4614, 8962, 1791, 1818,
                                                                       4884, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9472, 0, 3, 6,
                                                                       7292, 3630, 7322, 8962,
                                                                       4632, 8992, 1818, 1845,
                                                                       4938, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 9562, 0, 3, 6,
                                                                       7352, 3666, 7412, 9022,
                                                                       4668, 9112, 1899, 1953,
                                                                       4992, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 9742, 0, 3, 6,
                                                                       7412, 3702, 7472, 9112,
                                                                       4722, 9202, 1953, 2007,
                                                                       5100, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 9922, 0, 3, 6,
                                                                       7472, 3738, 7532, 9202,
                                                                       4776, 9292, 2007, 2061,
                                                                       5208, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 10102, 0, 3, 6,
                                                                       7532, 3774, 7592, 9292,
                                                                       4830, 9382, 2061, 2115,
                                                                       5316, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 10282, 0, 3, 6,
                                                                       7592, 3810, 7652, 9382,
                                                                       4884, 9472, 2115, 2169,
                                                                       5424, ncols, gamma, p,
                                                                       q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 10462, 0, 3, 6,
                                                                       7712, 3882, 7812, 9022,
                                                                       9112, 9562, 4992, 9742,
                                                                       2277, 2367, 5532, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 10762, 0, 3, 6,
                                                                       7812, 3942, 7912, 9112,
                                                                       9202, 9742, 5100, 9922,
                                                                       2367, 2457, 5712, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 11062, 0, 3, 6,
                                                                       7912, 4002, 8012, 9202,
                                                                       9292, 9922, 5208, 10102,
                                                                       2457, 2547, 5892, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 11362, 0, 3, 6,
                                                                       8012, 4062, 8112, 9292,
                                                                       9382, 10102, 5316, 10282,
                                                                       2547, 2637, 6072, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 11662, 0, 3, 6,
                                                                       8212, 4182, 8362, 9562,
                                                                       9742, 10462, 5532, 10762,
                                                                       2817, 2952, 6252, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 12112, 0, 3, 6,
                                                                       8362, 4272, 8512, 9742,
                                                                       9922, 10762, 5712, 11062,
                                                                       2952, 3087, 6522, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 12562, 0, 3, 6,
                                                                       8512, 4362, 8662, 9922,
                                                                       10102, 11062, 5892, 11362,
                                                                       3087, 3222, 6792, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13012, 6, 3492,
                                                                       3498, 7082, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13027, 6, 3498,
                                                                       3504, 7092, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13042, 6, 3504,
                                                                       3510, 7102, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13057, 6, 3510,
                                                                       3516, 7112, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13072, 6, 3516,
                                                                       3522, 7122, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13087, 6, 3522,
                                                                       3528, 7132, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13102, 3, 6,
                                                                       13012, 7082, 13027, 3540,
                                                                       3558, 7202, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13147, 3, 6,
                                                                       13027, 7092, 13042, 3558,
                                                                       3576, 7232, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13192, 3, 6,
                                                                       13042, 7102, 13057, 3576,
                                                                       3594, 7262, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13237, 3, 6,
                                                                       13057, 7112, 13072, 3594,
                                                                       3612, 7292, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13282, 3, 6,
                                                                       13072, 7122, 13087, 3612,
                                                                       3630, 7322, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 13327, 3, 6,
                                                                       13102, 7202, 13147, 3666,
                                                                       3702, 7472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 13417, 3, 6,
                                                                       13147, 7232, 13192, 3702,
                                                                       3738, 7532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 13507, 3, 6,
                                                                       13192, 7262, 13237, 3738,
                                                                       3774, 7592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 13597, 3, 6,
                                                                       13237, 7292, 13282, 3774,
                                                                       3810, 7652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 13687, 3, 6,
                                                                       13327, 7472, 13417, 3882,
                                                                       3942, 7912, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 13837, 3, 6,
                                                                       13417, 7532, 13507, 3942,
                                                                       4002, 8012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 13987, 3, 6,
                                                                       13507, 7592, 13597, 4002,
                                                                       4062, 8112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 14137, 3, 6,
                                                                       13687, 7912, 13837, 4182,
                                                                       4272, 8512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 14362, 3, 6,
                                                                       13837, 8012, 13987, 4272,
                                                                       4362, 8662, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14587, 0, 6,
                                                                       13012, 7082, 13027, 4542,
                                                                       4560, 8872, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14632, 0, 6,
                                                                       13027, 7092, 13042, 4560,
                                                                       4578, 8902, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14677, 0, 6,
                                                                       13042, 7102, 13057, 4578,
                                                                       4596, 8932, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14722, 0, 6,
                                                                       13057, 7112, 13072, 4596,
                                                                       4614, 8962, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14767, 0, 6,
                                                                       13072, 7122, 13087, 4614,
                                                                       4632, 8992, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 14812, 0, 3, 6,
                                                                       13102, 7202, 13147, 14587,
                                                                       8872, 14632, 4668, 4722,
                                                                       9202, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 14947, 0, 3, 6,
                                                                       13147, 7232, 13192, 14632,
                                                                       8902, 14677, 4722, 4776,
                                                                       9292, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 15082, 0, 3, 6,
                                                                       13192, 7262, 13237, 14677,
                                                                       8932, 14722, 4776, 4830,
                                                                       9382, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 15217, 0, 3, 6,
                                                                       13237, 7292, 13282, 14722,
                                                                       8962, 14767, 4830, 4884,
                                                                       9472, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 15352, 0, 3, 6,
                                                                       13327, 7472, 13417, 14812,
                                                                       9202, 14947, 4992, 5100,
                                                                       9922, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 15622, 0, 3, 6,
                                                                       13417, 7532, 13507, 14947,
                                                                       9292, 15082, 5100, 5208,
                                                                       10102, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 15892, 0, 3, 6,
                                                                       13507, 7592, 13597, 15082,
                                                                       9382, 15217, 5208, 5316,
                                                                       10282, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 16162, 0, 3, 6,
                                                                       13687, 7912, 13837, 14812,
                                                                       14947, 15352, 9922, 15622,
                                                                       5532, 5712, 11062, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 16612, 0, 3, 6,
                                                                       13837, 8012, 13987, 14947,
                                                                       15082, 15622, 10102,
                                                                       15892, 5712, 5892, 11362,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 17062, 0, 3, 6,
                                                                       14137, 8512, 14362, 15352,
                                                                       15622, 16162, 11062,
                                                                       16612, 6252, 6522, 12562,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17737, 6, 7062,
                                                                       7072, 13012, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17758, 6, 7072,
                                                                       7082, 13027, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17779, 6, 7082,
                                                                       7092, 13042, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17800, 6, 7092,
                                                                       7102, 13057, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17821, 6, 7102,
                                                                       7112, 13072, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17842, 6, 7112,
                                                                       7122, 13087, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17863, 3, 6,
                                                                       17737, 13012, 17758, 7142,
                                                                       7172, 13102, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17926, 3, 6,
                                                                       17758, 13027, 17779, 7172,
                                                                       7202, 13147, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17989, 3, 6,
                                                                       17779, 13042, 17800, 7202,
                                                                       7232, 13192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 18052, 3, 6,
                                                                       17800, 13057, 17821, 7232,
                                                                       7262, 13237, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 18115, 3, 6,
                                                                       17821, 13072, 17842, 7262,
                                                                       7292, 13282, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 18178, 3, 6,
                                                                       17863, 13102, 17926, 7352,
                                                                       7412, 13327, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 18304, 3, 6,
                                                                       17926, 13147, 17989, 7412,
                                                                       7472, 13417, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 18430, 3, 6,
                                                                       17989, 13192, 18052, 7472,
                                                                       7532, 13507, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 18556, 3, 6,
                                                                       18052, 13237, 18115, 7532,
                                                                       7592, 13597, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 18682, 3, 6,
                                                                       18178, 13327, 18304, 7712,
                                                                       7812, 13687, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 18892, 3, 6,
                                                                       18304, 13417, 18430, 7812,
                                                                       7912, 13837, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 19102, 3, 6,
                                                                       18430, 13507, 18556, 7912,
                                                                       8012, 13987, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 19312, 3, 6,
                                                                       18682, 13687, 18892, 8212,
                                                                       8362, 14137, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 19627, 3, 6,
                                                                       18892, 13837, 19102, 8362,
                                                                       8512, 14362, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 19942, 0, 6,
                                                                       17737, 13012, 17758, 8812,
                                                                       8842, 14587, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 20005, 0, 6,
                                                                       17758, 13027, 17779, 8842,
                                                                       8872, 14632, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 20068, 0, 6,
                                                                       17779, 13042, 17800, 8872,
                                                                       8902, 14677, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 20131, 0, 6,
                                                                       17800, 13057, 17821, 8902,
                                                                       8932, 14722, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 20194, 0, 6,
                                                                       17821, 13072, 17842, 8932,
                                                                       8962, 14767, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 20257, 0, 3, 6,
                                                                       17863, 13102, 17926,
                                                                       19942, 14587, 20005, 9022,
                                                                       9112, 14812, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 20446, 0, 3, 6,
                                                                       17926, 13147, 17989,
                                                                       20005, 14632, 20068, 9112,
                                                                       9202, 14947, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 20635, 0, 3, 6,
                                                                       17989, 13192, 18052,
                                                                       20068, 14677, 20131, 9202,
                                                                       9292, 15082, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 20824, 0, 3, 6,
                                                                       18052, 13237, 18115,
                                                                       20131, 14722, 20194, 9292,
                                                                       9382, 15217, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 21013, 0, 3, 6,
                                                                       18178, 13327, 18304,
                                                                       20257, 14812, 20446, 9562,
                                                                       9742, 15352, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 21391, 0, 3, 6,
                                                                       18304, 13417, 18430,
                                                                       20446, 14947, 20635, 9742,
                                                                       9922, 15622, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 21769, 0, 3, 6,
                                                                       18430, 13507, 18556,
                                                                       20635, 15082, 20824, 9922,
                                                                       10102, 15892, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 22147, 0, 3, 6,
                                                                       18682, 13687, 18892,
                                                                       20257, 20446, 21013,
                                                                       15352, 21391, 10462,
                                                                       10762, 16162, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 22777, 0, 3, 6,
                                                                       18892, 13837, 19102,
                                                                       20446, 20635, 21391,
                                                                       15622, 21769, 10762,
                                                                       11062, 16612, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgh_three_center_electron_repulsion_0(buffer, 23407, 0, 3, 6,
                                                                       19312, 14137, 19627,
                                                                       21013, 21391, 22147,
                                                                       16162, 22777, 11662,
                                                                       12112, 17062, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 24352, 23407, 1, 315, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 24667, 23407, 1, 315, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 24982, 23407, 1, 315, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 25297, 24352, 945, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 26242, 25297, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 26242, 11, nmax);

        simdtrf::transform_h_inner(buffer, 26242, 25612, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 99 * nvalues + n * npairs, nvalues, buffer, 26242,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 26242, 25927, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 198 * nvalues + n * npairs, nvalues, buffer, 26242,
                                   11, nmax);
    }

    for (size_t m = 0; m < 297; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
