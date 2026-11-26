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


#include "SimdThreeCenterElectronRepulsionGeom100RecSGI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPDI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_sgi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_sgi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 42059, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 351 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 42059, 40604, 1260, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 9, 6, 11,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 3, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 3, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 55, 3, 6, 10, 11,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 61, 3, 6, 11, 12,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 67, 3, 6, 12, 13,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 73, 3, 6, 13, 14,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 79, 3, 6, 14, 15,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 85, 3, 6, 15, 16,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 91, 3, 6, 16, 17,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 97, 3, 6, 17, 18,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 103, 3, 6, 18, 19,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 109, 3, 6, 19, 20,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 115, 3, 6, 22, 25,
                                                                       55, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 125, 3, 6, 25, 28,
                                                                       61, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 135, 3, 6, 28, 31,
                                                                       67, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 145, 3, 6, 31, 34,
                                                                       73, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 155, 3, 6, 34, 37,
                                                                       79, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 165, 3, 6, 37, 40,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 175, 3, 6, 40, 43,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 185, 3, 6, 43, 46,
                                                                       97, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 195, 3, 6, 46, 49,
                                                                       103, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 205, 3, 6, 55, 61,
                                                                       115, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 220, 3, 6, 61, 67,
                                                                       125, 135, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 235, 3, 6, 67, 73,
                                                                       135, 145, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 250, 3, 6, 73, 79,
                                                                       145, 155, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 265, 3, 6, 79, 85,
                                                                       155, 165, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 280, 3, 6, 85, 91,
                                                                       165, 175, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 295, 3, 6, 91, 97,
                                                                       175, 185, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 310, 3, 6, 97,
                                                                       103, 185, 195, ncols,
                                                                       gamma, p, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 325, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 328, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 331, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 334, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 337, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 340, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 343, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 346, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 349, 0, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 352, 0, 6, 10, 11,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 361, 0, 6, 11, 12,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 370, 0, 6, 12, 13,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 379, 0, 6, 13, 14,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 388, 0, 6, 14, 15,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 397, 0, 6, 15, 16,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 406, 0, 6, 16, 17,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 415, 0, 6, 17, 18,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 424, 0, 6, 18, 19,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 433, 0, 6, 19, 20,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 442, 0, 3, 6, 22,
                                                                       25, 55, 61, 352, 361,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 460, 0, 3, 6, 25,
                                                                       28, 61, 67, 361, 370,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 478, 0, 3, 6, 28,
                                                                       31, 67, 73, 370, 379,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 496, 0, 3, 6, 31,
                                                                       34, 73, 79, 379, 388,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 514, 0, 3, 6, 34,
                                                                       37, 79, 85, 388, 397,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 532, 0, 3, 6, 37,
                                                                       40, 85, 91, 397, 406,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 550, 0, 3, 6, 40,
                                                                       43, 91, 97, 406, 415,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 568, 0, 3, 6, 43,
                                                                       46, 97, 103, 415, 424,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 586, 0, 3, 6, 46,
                                                                       49, 103, 109, 424, 433,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 604, 0, 3, 6, 55,
                                                                       61, 115, 125, 442, 460,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 634, 0, 3, 6, 61,
                                                                       67, 125, 135, 460, 478,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 664, 0, 3, 6, 67,
                                                                       73, 135, 145, 478, 496,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 694, 0, 3, 6, 73,
                                                                       79, 145, 155, 496, 514,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 724, 0, 3, 6, 79,
                                                                       85, 155, 165, 514, 532,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 754, 0, 3, 6, 85,
                                                                       91, 165, 175, 532, 550,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 784, 0, 3, 6, 91,
                                                                       97, 175, 185, 550, 568,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 814, 0, 3, 6, 97,
                                                                       103, 185, 195, 568, 586,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 844, 0, 3, 6, 115,
                                                                       125, 205, 220, 604, 634,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 889, 0, 3, 6, 125,
                                                                       135, 220, 235, 634, 664,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 934, 0, 3, 6, 135,
                                                                       145, 235, 250, 664, 694,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 979, 0, 3, 6, 145,
                                                                       155, 250, 265, 694, 724,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1024, 0, 3, 6,
                                                                       155, 165, 265, 280, 724,
                                                                       754, ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1069, 0, 3, 6,
                                                                       165, 175, 280, 295, 754,
                                                                       784, ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1114, 0, 3, 6,
                                                                       175, 185, 295, 310, 784,
                                                                       814, ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1159, 6, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1162, 6, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1165, 6, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1168, 6, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1171, 6, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1174, 6, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1177, 6, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1180, 6, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1183, 6, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1186, 6, 21,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1189, 6, 12, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1198, 6, 13, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1207, 6, 14, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1216, 6, 15, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1225, 6, 16, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1234, 6, 17, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1243, 6, 18, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1252, 6, 19, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1261, 6, 20, 52,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1270, 6, 28, 67,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1288, 6, 31, 73,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1306, 6, 34, 79,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1324, 6, 37, 85,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1342, 6, 40, 91,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1360, 6, 43, 97,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1378, 6, 46, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1396, 6, 49, 109,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1414, 6, 67, 135,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1444, 6, 73, 145,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1474, 6, 79, 155,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1504, 6, 85, 165,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1534, 6, 91, 175,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1564, 6, 97, 185,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1594, 6, 103, 195,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1624, 6, 135, 235,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1669, 6, 145, 250,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1714, 6, 155, 265,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1759, 6, 165, 280,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1804, 6, 175, 295,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1849, 6, 185, 310,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1894, 6, 12, 325,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1903, 6, 13, 328,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1912, 6, 14, 331,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1921, 6, 15, 334,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1930, 6, 16, 337,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1939, 6, 17, 340,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1948, 6, 18, 343,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1957, 6, 19, 346,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1966, 6, 20, 349,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1975, 6, 28, 325,
                                                                       370, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2002, 6, 31, 328,
                                                                       379, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2029, 6, 34, 331,
                                                                       388, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2056, 6, 37, 334,
                                                                       397, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2083, 6, 40, 337,
                                                                       406, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2110, 6, 43, 340,
                                                                       415, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2137, 6, 46, 343,
                                                                       424, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2164, 6, 49, 346,
                                                                       433, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2191, 3, 6, 67,
                                                                       1975, 370, 2002, 478,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2245, 3, 6, 73,
                                                                       2002, 379, 2029, 496,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2299, 3, 6, 79,
                                                                       2029, 388, 2056, 514,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2353, 3, 6, 85,
                                                                       2056, 397, 2083, 532,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2407, 3, 6, 91,
                                                                       2083, 406, 2110, 550,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2461, 3, 6, 97,
                                                                       2110, 415, 2137, 568,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2515, 3, 6, 103,
                                                                       2137, 424, 2164, 586,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2569, 3, 6, 135,
                                                                       2191, 478, 2245, 664,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2659, 3, 6, 145,
                                                                       2245, 496, 2299, 694,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2749, 3, 6, 155,
                                                                       2299, 514, 2353, 724,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2839, 3, 6, 165,
                                                                       2353, 532, 2407, 754,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2929, 3, 6, 175,
                                                                       2407, 550, 2461, 784,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3019, 3, 6, 185,
                                                                       2461, 568, 2515, 814,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 3109, 3, 6, 235,
                                                                       2569, 664, 2659, 934,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 3244, 3, 6, 250,
                                                                       2659, 694, 2749, 979,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 3379, 3, 6, 265,
                                                                       2749, 724, 2839, 1024,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 3514, 3, 6, 280,
                                                                       2839, 754, 2929, 1069,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 3649, 3, 6, 295,
                                                                       2929, 784, 3019, 1114,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3784, 6, 10, 11,
                                                                       1159, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3790, 6, 11, 12,
                                                                       1162, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3796, 6, 12, 13,
                                                                       1165, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3802, 6, 13, 14,
                                                                       1168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3808, 6, 14, 15,
                                                                       1171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3814, 6, 15, 16,
                                                                       1174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3820, 6, 16, 17,
                                                                       1177, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3826, 6, 17, 18,
                                                                       1180, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3832, 6, 18, 19,
                                                                       1183, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3838, 6, 19, 20,
                                                                       1186, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3844, 3, 6, 3784,
                                                                       1159, 3790, 1189, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3862, 3, 6, 3790,
                                                                       1162, 3796, 1198, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3880, 3, 6, 3796,
                                                                       1165, 3802, 1207, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3898, 3, 6, 3802,
                                                                       1168, 3808, 1216, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3916, 3, 6, 3808,
                                                                       1171, 3814, 1225, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3934, 3, 6, 3814,
                                                                       1174, 3820, 1234, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3952, 3, 6, 3820,
                                                                       1177, 3826, 1243, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3970, 3, 6, 3826,
                                                                       1180, 3832, 1252, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3988, 3, 6, 3832,
                                                                       1183, 3838, 1261, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4006, 3, 6, 3844,
                                                                       1189, 3862, 55, 61, 1270,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4042, 3, 6, 3862,
                                                                       1198, 3880, 61, 67, 1288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4078, 3, 6, 3880,
                                                                       1207, 3898, 67, 73, 1306,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4114, 3, 6, 3898,
                                                                       1216, 3916, 73, 79, 1324,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4150, 3, 6, 3916,
                                                                       1225, 3934, 79, 85, 1342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4186, 3, 6, 3934,
                                                                       1234, 3952, 85, 91, 1360,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4222, 3, 6, 3952,
                                                                       1243, 3970, 91, 97, 1378,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4258, 3, 6, 3970,
                                                                       1252, 3988, 97, 103, 1396,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4294, 3, 6, 4006,
                                                                       1270, 4042, 115, 125,
                                                                       1414, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4354, 3, 6, 4042,
                                                                       1288, 4078, 125, 135,
                                                                       1444, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4414, 3, 6, 4078,
                                                                       1306, 4114, 135, 145,
                                                                       1474, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4474, 3, 6, 4114,
                                                                       1324, 4150, 145, 155,
                                                                       1504, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4534, 3, 6, 4150,
                                                                       1342, 4186, 155, 165,
                                                                       1534, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4594, 3, 6, 4186,
                                                                       1360, 4222, 165, 175,
                                                                       1564, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4654, 3, 6, 4222,
                                                                       1378, 4258, 175, 185,
                                                                       1594, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4714, 3, 6, 4294,
                                                                       1414, 4354, 205, 220,
                                                                       1624, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4804, 3, 6, 4354,
                                                                       1444, 4414, 220, 235,
                                                                       1669, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4894, 3, 6, 4414,
                                                                       1474, 4474, 235, 250,
                                                                       1714, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 4984, 3, 6, 4474,
                                                                       1504, 4534, 250, 265,
                                                                       1759, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5074, 3, 6, 4534,
                                                                       1534, 4594, 265, 280,
                                                                       1804, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5164, 3, 6, 4594,
                                                                       1564, 4654, 280, 295,
                                                                       1849, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5254, 0, 6, 3784,
                                                                       1159, 3790, 1894, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5272, 0, 6, 3790,
                                                                       1162, 3796, 1903, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5290, 0, 6, 3796,
                                                                       1165, 3802, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5308, 0, 6, 3802,
                                                                       1168, 3808, 1921, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5326, 0, 6, 3808,
                                                                       1171, 3814, 1930, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5344, 0, 6, 3814,
                                                                       1174, 3820, 1939, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5362, 0, 6, 3820,
                                                                       1177, 3826, 1948, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5380, 0, 6, 3826,
                                                                       1180, 3832, 1957, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5398, 0, 6, 3832,
                                                                       1183, 3838, 1966, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5416, 0, 3, 6,
                                                                       3844, 1189, 3862, 5254,
                                                                       1894, 5272, 352, 361,
                                                                       1975, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5470, 0, 3, 6,
                                                                       3862, 1198, 3880, 5272,
                                                                       1903, 5290, 361, 370,
                                                                       2002, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5524, 0, 3, 6,
                                                                       3880, 1207, 3898, 5290,
                                                                       1912, 5308, 370, 379,
                                                                       2029, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5578, 0, 3, 6,
                                                                       3898, 1216, 3916, 5308,
                                                                       1921, 5326, 379, 388,
                                                                       2056, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5632, 0, 3, 6,
                                                                       3916, 1225, 3934, 5326,
                                                                       1930, 5344, 388, 397,
                                                                       2083, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5686, 0, 3, 6,
                                                                       3934, 1234, 3952, 5344,
                                                                       1939, 5362, 397, 406,
                                                                       2110, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5740, 0, 3, 6,
                                                                       3952, 1243, 3970, 5362,
                                                                       1948, 5380, 406, 415,
                                                                       2137, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5794, 0, 3, 6,
                                                                       3970, 1252, 3988, 5380,
                                                                       1957, 5398, 415, 424,
                                                                       2164, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5848, 0, 3, 6,
                                                                       4006, 1270, 4042, 5416,
                                                                       1975, 5470, 442, 460,
                                                                       2191, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5956, 0, 3, 6,
                                                                       4042, 1288, 4078, 5470,
                                                                       2002, 5524, 460, 478,
                                                                       2245, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6064, 0, 3, 6,
                                                                       4078, 1306, 4114, 5524,
                                                                       2029, 5578, 478, 496,
                                                                       2299, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6172, 0, 3, 6,
                                                                       4114, 1324, 4150, 5578,
                                                                       2056, 5632, 496, 514,
                                                                       2353, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6280, 0, 3, 6,
                                                                       4150, 1342, 4186, 5632,
                                                                       2083, 5686, 514, 532,
                                                                       2407, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6388, 0, 3, 6,
                                                                       4186, 1360, 4222, 5686,
                                                                       2110, 5740, 532, 550,
                                                                       2461, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6496, 0, 3, 6,
                                                                       4222, 1378, 4258, 5740,
                                                                       2137, 5794, 550, 568,
                                                                       2515, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 6604, 0, 3, 6,
                                                                       4294, 1414, 4354, 5416,
                                                                       5470, 5848, 2191, 5956,
                                                                       604, 634, 2569, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 6784, 0, 3, 6,
                                                                       4354, 1444, 4414, 5470,
                                                                       5524, 5956, 2245, 6064,
                                                                       634, 664, 2659, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 6964, 0, 3, 6,
                                                                       4414, 1474, 4474, 5524,
                                                                       5578, 6064, 2299, 6172,
                                                                       664, 694, 2749, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 7144, 0, 3, 6,
                                                                       4474, 1504, 4534, 5578,
                                                                       5632, 6172, 2353, 6280,
                                                                       694, 724, 2839, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 7324, 0, 3, 6,
                                                                       4534, 1534, 4594, 5632,
                                                                       5686, 6280, 2407, 6388,
                                                                       724, 754, 2929, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 7504, 0, 3, 6,
                                                                       4594, 1564, 4654, 5686,
                                                                       5740, 6388, 2461, 6496,
                                                                       754, 784, 3019, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 7684, 0, 3, 6,
                                                                       4714, 1624, 4804, 5848,
                                                                       5956, 6604, 2569, 6784,
                                                                       844, 889, 3109, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 7954, 0, 3, 6,
                                                                       4804, 1669, 4894, 5956,
                                                                       6064, 6784, 2659, 6964,
                                                                       889, 934, 3244, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 8224, 0, 3, 6,
                                                                       4894, 1714, 4984, 6064,
                                                                       6172, 6964, 2749, 7144,
                                                                       934, 979, 3379, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 8494, 0, 3, 6,
                                                                       4984, 1759, 5074, 6172,
                                                                       6280, 7144, 2839, 7324,
                                                                       979, 1024, 3514, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 8764, 0, 3, 6,
                                                                       5074, 1804, 5164, 6280,
                                                                       6388, 7324, 2929, 7504,
                                                                       1024, 1069, 3649, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9034, 6, 1159,
                                                                       1162, 3796, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9044, 6, 1162,
                                                                       1165, 3802, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9054, 6, 1165,
                                                                       1168, 3808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9064, 6, 1168,
                                                                       1171, 3814, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9074, 6, 1171,
                                                                       1174, 3820, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9084, 6, 1174,
                                                                       1177, 3826, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9094, 6, 1177,
                                                                       1180, 3832, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9104, 6, 1180,
                                                                       1183, 3838, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9114, 3, 6, 9034,
                                                                       3796, 9044, 3880, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9144, 3, 6, 9044,
                                                                       3802, 9054, 3898, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9174, 3, 6, 9054,
                                                                       3808, 9064, 3916, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9204, 3, 6, 9064,
                                                                       3814, 9074, 3934, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9234, 3, 6, 9074,
                                                                       3820, 9084, 3952, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9264, 3, 6, 9084,
                                                                       3826, 9094, 3970, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9294, 3, 6, 9094,
                                                                       3832, 9104, 3988, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9324, 3, 6, 9114,
                                                                       3880, 9144, 1270, 1288,
                                                                       4078, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9384, 3, 6, 9144,
                                                                       3898, 9174, 1288, 1306,
                                                                       4114, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9444, 3, 6, 9174,
                                                                       3916, 9204, 1306, 1324,
                                                                       4150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9504, 3, 6, 9204,
                                                                       3934, 9234, 1324, 1342,
                                                                       4186, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9564, 3, 6, 9234,
                                                                       3952, 9264, 1342, 1360,
                                                                       4222, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9624, 3, 6, 9264,
                                                                       3970, 9294, 1360, 1378,
                                                                       4258, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9684, 3, 6, 9324,
                                                                       4078, 9384, 1414, 1444,
                                                                       4414, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9784, 3, 6, 9384,
                                                                       4114, 9444, 1444, 1474,
                                                                       4474, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9884, 3, 6, 9444,
                                                                       4150, 9504, 1474, 1504,
                                                                       4534, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9984, 3, 6, 9504,
                                                                       4186, 9564, 1504, 1534,
                                                                       4594, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10084, 3, 6, 9564,
                                                                       4222, 9624, 1534, 1564,
                                                                       4654, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10184, 3, 6, 9684,
                                                                       4414, 9784, 1624, 1669,
                                                                       4894, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10334, 3, 6, 9784,
                                                                       4474, 9884, 1669, 1714,
                                                                       4984, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10484, 3, 6, 9884,
                                                                       4534, 9984, 1714, 1759,
                                                                       5074, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10634, 3, 6, 9984,
                                                                       4594, 10084, 1759, 1804,
                                                                       5164, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10784, 0, 6, 9034,
                                                                       3796, 9044, 1894, 1903,
                                                                       5290, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10814, 0, 6, 9044,
                                                                       3802, 9054, 1903, 1912,
                                                                       5308, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10844, 0, 6, 9054,
                                                                       3808, 9064, 1912, 1921,
                                                                       5326, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10874, 0, 6, 9064,
                                                                       3814, 9074, 1921, 1930,
                                                                       5344, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10904, 0, 6, 9074,
                                                                       3820, 9084, 1930, 1939,
                                                                       5362, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10934, 0, 6, 9084,
                                                                       3826, 9094, 1939, 1948,
                                                                       5380, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10964, 0, 6, 9094,
                                                                       3832, 9104, 1948, 1957,
                                                                       5398, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10994, 0, 3, 6,
                                                                       9114, 3880, 9144, 10784,
                                                                       5290, 10814, 1975, 2002,
                                                                       5524, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11084, 0, 3, 6,
                                                                       9144, 3898, 9174, 10814,
                                                                       5308, 10844, 2002, 2029,
                                                                       5578, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11174, 0, 3, 6,
                                                                       9174, 3916, 9204, 10844,
                                                                       5326, 10874, 2029, 2056,
                                                                       5632, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11264, 0, 3, 6,
                                                                       9204, 3934, 9234, 10874,
                                                                       5344, 10904, 2056, 2083,
                                                                       5686, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11354, 0, 3, 6,
                                                                       9234, 3952, 9264, 10904,
                                                                       5362, 10934, 2083, 2110,
                                                                       5740, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11444, 0, 3, 6,
                                                                       9264, 3970, 9294, 10934,
                                                                       5380, 10964, 2110, 2137,
                                                                       5794, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 11534, 0, 3, 6,
                                                                       9324, 4078, 9384, 10994,
                                                                       5524, 11084, 2191, 2245,
                                                                       6064, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 11714, 0, 3, 6,
                                                                       9384, 4114, 9444, 11084,
                                                                       5578, 11174, 2245, 2299,
                                                                       6172, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 11894, 0, 3, 6,
                                                                       9444, 4150, 9504, 11174,
                                                                       5632, 11264, 2299, 2353,
                                                                       6280, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 12074, 0, 3, 6,
                                                                       9504, 4186, 9564, 11264,
                                                                       5686, 11354, 2353, 2407,
                                                                       6388, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 12254, 0, 3, 6,
                                                                       9564, 4222, 9624, 11354,
                                                                       5740, 11444, 2407, 2461,
                                                                       6496, ncols, gamma, p,
                                                                       q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 12434, 0, 3, 6,
                                                                       9684, 4414, 9784, 10994,
                                                                       11084, 11534, 6064, 11714,
                                                                       2569, 2659, 6964, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 12734, 0, 3, 6,
                                                                       9784, 4474, 9884, 11084,
                                                                       11174, 11714, 6172, 11894,
                                                                       2659, 2749, 7144, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 13034, 0, 3, 6,
                                                                       9884, 4534, 9984, 11174,
                                                                       11264, 11894, 6280, 12074,
                                                                       2749, 2839, 7324, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 13334, 0, 3, 6,
                                                                       9984, 4594, 10084, 11264,
                                                                       11354, 12074, 6388, 12254,
                                                                       2839, 2929, 7504, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 13634, 0, 3, 6,
                                                                       10184, 4894, 10334, 11534,
                                                                       11714, 12434, 6964, 12734,
                                                                       3109, 3244, 8224, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 14084, 0, 3, 6,
                                                                       10334, 4984, 10484, 11714,
                                                                       11894, 12734, 7144, 13034,
                                                                       3244, 3379, 8494, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 14534, 0, 3, 6,
                                                                       10484, 5074, 10634, 11894,
                                                                       12074, 13034, 7324, 13334,
                                                                       3379, 3514, 8764, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14984, 6, 3784,
                                                                       3790, 9034, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14999, 6, 3790,
                                                                       3796, 9044, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15014, 6, 3796,
                                                                       3802, 9054, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15029, 6, 3802,
                                                                       3808, 9064, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15044, 6, 3808,
                                                                       3814, 9074, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15059, 6, 3814,
                                                                       3820, 9084, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15074, 6, 3820,
                                                                       3826, 9094, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15089, 6, 3826,
                                                                       3832, 9104, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15104, 3, 6,
                                                                       14984, 9034, 14999, 3844,
                                                                       3862, 9114, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15149, 3, 6,
                                                                       14999, 9044, 15014, 3862,
                                                                       3880, 9144, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15194, 3, 6,
                                                                       15014, 9054, 15029, 3880,
                                                                       3898, 9174, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15239, 3, 6,
                                                                       15029, 9064, 15044, 3898,
                                                                       3916, 9204, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15284, 3, 6,
                                                                       15044, 9074, 15059, 3916,
                                                                       3934, 9234, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15329, 3, 6,
                                                                       15059, 9084, 15074, 3934,
                                                                       3952, 9264, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15374, 3, 6,
                                                                       15074, 9094, 15089, 3952,
                                                                       3970, 9294, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 15419, 3, 6,
                                                                       15104, 9114, 15149, 4006,
                                                                       4042, 9324, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 15509, 3, 6,
                                                                       15149, 9144, 15194, 4042,
                                                                       4078, 9384, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 15599, 3, 6,
                                                                       15194, 9174, 15239, 4078,
                                                                       4114, 9444, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 15689, 3, 6,
                                                                       15239, 9204, 15284, 4114,
                                                                       4150, 9504, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 15779, 3, 6,
                                                                       15284, 9234, 15329, 4150,
                                                                       4186, 9564, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 15869, 3, 6,
                                                                       15329, 9264, 15374, 4186,
                                                                       4222, 9624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 15959, 3, 6,
                                                                       15419, 9324, 15509, 4294,
                                                                       4354, 9684, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 16109, 3, 6,
                                                                       15509, 9384, 15599, 4354,
                                                                       4414, 9784, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 16259, 3, 6,
                                                                       15599, 9444, 15689, 4414,
                                                                       4474, 9884, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 16409, 3, 6,
                                                                       15689, 9504, 15779, 4474,
                                                                       4534, 9984, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 16559, 3, 6,
                                                                       15779, 9564, 15869, 4534,
                                                                       4594, 10084, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 16709, 3, 6,
                                                                       15959, 9684, 16109, 4714,
                                                                       4804, 10184, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 16934, 3, 6,
                                                                       16109, 9784, 16259, 4804,
                                                                       4894, 10334, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 17159, 3, 6,
                                                                       16259, 9884, 16409, 4894,
                                                                       4984, 10484, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 17384, 3, 6,
                                                                       16409, 9984, 16559, 4984,
                                                                       5074, 10634, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17609, 0, 6,
                                                                       14984, 9034, 14999, 5254,
                                                                       5272, 10784, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17654, 0, 6,
                                                                       14999, 9044, 15014, 5272,
                                                                       5290, 10814, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17699, 0, 6,
                                                                       15014, 9054, 15029, 5290,
                                                                       5308, 10844, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17744, 0, 6,
                                                                       15029, 9064, 15044, 5308,
                                                                       5326, 10874, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17789, 0, 6,
                                                                       15044, 9074, 15059, 5326,
                                                                       5344, 10904, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17834, 0, 6,
                                                                       15059, 9084, 15074, 5344,
                                                                       5362, 10934, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17879, 0, 6,
                                                                       15074, 9094, 15089, 5362,
                                                                       5380, 10964, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 17924, 0, 3, 6,
                                                                       15104, 9114, 15149, 17609,
                                                                       10784, 17654, 5416, 5470,
                                                                       10994, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 18059, 0, 3, 6,
                                                                       15149, 9144, 15194, 17654,
                                                                       10814, 17699, 5470, 5524,
                                                                       11084, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 18194, 0, 3, 6,
                                                                       15194, 9174, 15239, 17699,
                                                                       10844, 17744, 5524, 5578,
                                                                       11174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 18329, 0, 3, 6,
                                                                       15239, 9204, 15284, 17744,
                                                                       10874, 17789, 5578, 5632,
                                                                       11264, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 18464, 0, 3, 6,
                                                                       15284, 9234, 15329, 17789,
                                                                       10904, 17834, 5632, 5686,
                                                                       11354, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 18599, 0, 3, 6,
                                                                       15329, 9264, 15374, 17834,
                                                                       10934, 17879, 5686, 5740,
                                                                       11444, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 18734, 0, 3, 6,
                                                                       15419, 9324, 15509, 17924,
                                                                       10994, 18059, 5848, 5956,
                                                                       11534, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 19004, 0, 3, 6,
                                                                       15509, 9384, 15599, 18059,
                                                                       11084, 18194, 5956, 6064,
                                                                       11714, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 19274, 0, 3, 6,
                                                                       15599, 9444, 15689, 18194,
                                                                       11174, 18329, 6064, 6172,
                                                                       11894, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 19544, 0, 3, 6,
                                                                       15689, 9504, 15779, 18329,
                                                                       11264, 18464, 6172, 6280,
                                                                       12074, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 19814, 0, 3, 6,
                                                                       15779, 9564, 15869, 18464,
                                                                       11354, 18599, 6280, 6388,
                                                                       12254, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 20084, 0, 3, 6,
                                                                       15959, 9684, 16109, 17924,
                                                                       18059, 18734, 11534,
                                                                       19004, 6604, 6784, 12434,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 20534, 0, 3, 6,
                                                                       16109, 9784, 16259, 18059,
                                                                       18194, 19004, 11714,
                                                                       19274, 6784, 6964, 12734,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 20984, 0, 3, 6,
                                                                       16259, 9884, 16409, 18194,
                                                                       18329, 19274, 11894,
                                                                       19544, 6964, 7144, 13034,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 21434, 0, 3, 6,
                                                                       16409, 9984, 16559, 18329,
                                                                       18464, 19544, 12074,
                                                                       19814, 7144, 7324, 13334,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 21884, 0, 3, 6,
                                                                       16709, 10184, 16934,
                                                                       18734, 19004, 20084,
                                                                       12434, 20534, 7684, 7954,
                                                                       13634, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 22559, 0, 3, 6,
                                                                       16934, 10334, 17159,
                                                                       19004, 19274, 20534,
                                                                       12734, 20984, 7954, 8224,
                                                                       14084, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 23234, 0, 3, 6,
                                                                       17159, 10484, 17384,
                                                                       19274, 19544, 20984,
                                                                       13034, 21434, 8224, 8494,
                                                                       14534, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 23909, 6, 9034,
                                                                       9044, 15014, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 23930, 6, 9044,
                                                                       9054, 15029, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 23951, 6, 9054,
                                                                       9064, 15044, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 23972, 6, 9064,
                                                                       9074, 15059, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 23993, 6, 9074,
                                                                       9084, 15074, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 24014, 6, 9084,
                                                                       9094, 15089, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 24035, 3, 6,
                                                                       23909, 15014, 23930, 9114,
                                                                       9144, 15194, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 24098, 3, 6,
                                                                       23930, 15029, 23951, 9144,
                                                                       9174, 15239, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 24161, 3, 6,
                                                                       23951, 15044, 23972, 9174,
                                                                       9204, 15284, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 24224, 3, 6,
                                                                       23972, 15059, 23993, 9204,
                                                                       9234, 15329, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 24287, 3, 6,
                                                                       23993, 15074, 24014, 9234,
                                                                       9264, 15374, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 24350, 3, 6,
                                                                       24035, 15194, 24098, 9324,
                                                                       9384, 15599, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 24476, 3, 6,
                                                                       24098, 15239, 24161, 9384,
                                                                       9444, 15689, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 24602, 3, 6,
                                                                       24161, 15284, 24224, 9444,
                                                                       9504, 15779, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 24728, 3, 6,
                                                                       24224, 15329, 24287, 9504,
                                                                       9564, 15869, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 24854, 3, 6,
                                                                       24350, 15599, 24476, 9684,
                                                                       9784, 16259, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 25064, 3, 6,
                                                                       24476, 15689, 24602, 9784,
                                                                       9884, 16409, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 25274, 3, 6,
                                                                       24602, 15779, 24728, 9884,
                                                                       9984, 16559, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 25484, 3, 6,
                                                                       24854, 16259, 25064,
                                                                       10184, 10334, 17159,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 25799, 3, 6,
                                                                       25064, 16409, 25274,
                                                                       10334, 10484, 17384,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26114, 0, 6,
                                                                       23909, 15014, 23930,
                                                                       10784, 10814, 17699,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26177, 0, 6,
                                                                       23930, 15029, 23951,
                                                                       10814, 10844, 17744,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26240, 0, 6,
                                                                       23951, 15044, 23972,
                                                                       10844, 10874, 17789,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26303, 0, 6,
                                                                       23972, 15059, 23993,
                                                                       10874, 10904, 17834,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26366, 0, 6,
                                                                       23993, 15074, 24014,
                                                                       10904, 10934, 17879,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 26429, 0, 3, 6,
                                                                       24035, 15194, 24098,
                                                                       26114, 17699, 26177,
                                                                       10994, 11084, 18194,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 26618, 0, 3, 6,
                                                                       24098, 15239, 24161,
                                                                       26177, 17744, 26240,
                                                                       11084, 11174, 18329,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 26807, 0, 3, 6,
                                                                       24161, 15284, 24224,
                                                                       26240, 17789, 26303,
                                                                       11174, 11264, 18464,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 26996, 0, 3, 6,
                                                                       24224, 15329, 24287,
                                                                       26303, 17834, 26366,
                                                                       11264, 11354, 18599,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 27185, 0, 3, 6,
                                                                       24350, 15599, 24476,
                                                                       26429, 18194, 26618,
                                                                       11534, 11714, 19274,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 27563, 0, 3, 6,
                                                                       24476, 15689, 24602,
                                                                       26618, 18329, 26807,
                                                                       11714, 11894, 19544,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 27941, 0, 3, 6,
                                                                       24602, 15779, 24728,
                                                                       26807, 18464, 26996,
                                                                       11894, 12074, 19814,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 28319, 0, 3, 6,
                                                                       24854, 16259, 25064,
                                                                       26429, 26618, 27185,
                                                                       19274, 27563, 12434,
                                                                       12734, 20984, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 28949, 0, 3, 6,
                                                                       25064, 16409, 25274,
                                                                       26618, 26807, 27563,
                                                                       19544, 27941, 12734,
                                                                       13034, 21434, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgh_three_center_electron_repulsion_0(buffer, 29579, 0, 3, 6,
                                                                       25484, 17159, 25799,
                                                                       27185, 27563, 28319,
                                                                       20984, 28949, 13634,
                                                                       14084, 23234, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 30524, 6, 14984,
                                                                       14999, 23909, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 30552, 6, 14999,
                                                                       15014, 23930, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 30580, 6, 15014,
                                                                       15029, 23951, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 30608, 6, 15029,
                                                                       15044, 23972, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 30636, 6, 15044,
                                                                       15059, 23993, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 30664, 6, 15059,
                                                                       15074, 24014, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 30692, 3, 6,
                                                                       30524, 23909, 30552,
                                                                       15104, 15149, 24035,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 30776, 3, 6,
                                                                       30552, 23930, 30580,
                                                                       15149, 15194, 24098,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 30860, 3, 6,
                                                                       30580, 23951, 30608,
                                                                       15194, 15239, 24161,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 30944, 3, 6,
                                                                       30608, 23972, 30636,
                                                                       15239, 15284, 24224,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 31028, 3, 6,
                                                                       30636, 23993, 30664,
                                                                       15284, 15329, 24287,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 31112, 3, 6,
                                                                       30692, 24035, 30776,
                                                                       15419, 15509, 24350,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 31280, 3, 6,
                                                                       30776, 24098, 30860,
                                                                       15509, 15599, 24476,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 31448, 3, 6,
                                                                       30860, 24161, 30944,
                                                                       15599, 15689, 24602,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 31616, 3, 6,
                                                                       30944, 24224, 31028,
                                                                       15689, 15779, 24728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 31784, 3, 6,
                                                                       31112, 24350, 31280,
                                                                       15959, 16109, 24854,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 32064, 3, 6,
                                                                       31280, 24476, 31448,
                                                                       16109, 16259, 25064,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 32344, 3, 6,
                                                                       31448, 24602, 31616,
                                                                       16259, 16409, 25274,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 32624, 3, 6,
                                                                       31784, 24854, 32064,
                                                                       16709, 16934, 25484,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 33044, 3, 6,
                                                                       32064, 25064, 32344,
                                                                       16934, 17159, 25799,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 33464, 0, 6,
                                                                       30524, 23909, 30552,
                                                                       17609, 17654, 26114,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 33548, 0, 6,
                                                                       30552, 23930, 30580,
                                                                       17654, 17699, 26177,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 33632, 0, 6,
                                                                       30580, 23951, 30608,
                                                                       17699, 17744, 26240,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 33716, 0, 6,
                                                                       30608, 23972, 30636,
                                                                       17744, 17789, 26303,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 33800, 0, 6,
                                                                       30636, 23993, 30664,
                                                                       17789, 17834, 26366,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 33884, 0, 3, 6,
                                                                       30692, 24035, 30776,
                                                                       33464, 26114, 33548,
                                                                       17924, 18059, 26429,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 34136, 0, 3, 6,
                                                                       30776, 24098, 30860,
                                                                       33548, 26177, 33632,
                                                                       18059, 18194, 26618,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 34388, 0, 3, 6,
                                                                       30860, 24161, 30944,
                                                                       33632, 26240, 33716,
                                                                       18194, 18329, 26807,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 34640, 0, 3, 6,
                                                                       30944, 24224, 31028,
                                                                       33716, 26303, 33800,
                                                                       18329, 18464, 26996,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 34892, 0, 3, 6,
                                                                       31112, 24350, 31280,
                                                                       33884, 26429, 34136,
                                                                       18734, 19004, 27185,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 35396, 0, 3, 6,
                                                                       31280, 24476, 31448,
                                                                       34136, 26618, 34388,
                                                                       19004, 19274, 27563,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 35900, 0, 3, 6,
                                                                       31448, 24602, 31616,
                                                                       34388, 26807, 34640,
                                                                       19274, 19544, 27941,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfi_three_center_electron_repulsion_0(buffer, 36404, 0, 3, 6,
                                                                       31784, 24854, 32064,
                                                                       33884, 34136, 34892,
                                                                       27185, 35396, 20084,
                                                                       20534, 28319, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfi_three_center_electron_repulsion_0(buffer, 37244, 0, 3, 6,
                                                                       32064, 25064, 32344,
                                                                       34136, 34388, 35396,
                                                                       27563, 35900, 20534,
                                                                       20984, 28949, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgi_three_center_electron_repulsion_0(buffer, 38084, 0, 3, 6,
                                                                       32624, 25484, 33044,
                                                                       34892, 35396, 36404,
                                                                       28319, 37244, 21884,
                                                                       22559, 29579, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 39344, 38084, 1, 420, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 39764, 38084, 1, 420, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 40184, 38084, 1, 420, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 40604, 39344, 1260, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 41864, 40604, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 41864, 13, nmax);

        simdtrf::transform_i_inner(buffer, 41864, 41024, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 117 * nvalues + n * npairs, nvalues, buffer, 41864,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 41864, 41444, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 234 * nvalues + n * npairs, nvalues, buffer, 41864,
                                   13, nmax);
    }

    for (size_t m = 0; m < 351; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
