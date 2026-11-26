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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecSGI.hpp"

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
compute_rs_geom_100_sgi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_sgi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 83914, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 702 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 83914, 81199, 2520, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 9, 6, 11,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 22, 6, 11,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 56, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 59, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 62, 3, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 65, 3, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 68, 3, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 71, 3, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 74, 3, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 77, 3, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 80, 3, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 83, 3, 6, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 86, 3, 6, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 89, 3, 6, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 92, 3, 6, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 95, 3, 6, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 98, 3, 6, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 101, 3, 6, 10, 11,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 107, 3, 6, 11, 12,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 113, 3, 6, 12, 13,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 119, 3, 6, 13, 14,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 125, 3, 6, 14, 15,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 131, 3, 6, 15, 16,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 137, 3, 6, 16, 17,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 143, 3, 6, 17, 18,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 149, 3, 6, 18, 19,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 155, 3, 6, 19, 20,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 161, 3, 6, 23, 24,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 167, 3, 6, 24, 25,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 173, 3, 6, 25, 26,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 179, 3, 6, 26, 27,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 185, 3, 6, 27, 28,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 191, 3, 6, 28, 29,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 197, 3, 6, 29, 30,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 203, 3, 6, 30, 31,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 209, 3, 6, 31, 32,
                                                                       92, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 215, 3, 6, 32, 33,
                                                                       95, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 221, 3, 6, 35, 38,
                                                                       101, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 231, 3, 6, 38, 41,
                                                                       107, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 241, 3, 6, 41, 44,
                                                                       113, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 251, 3, 6, 44, 47,
                                                                       119, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 261, 3, 6, 47, 50,
                                                                       125, 131, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 271, 3, 6, 50, 53,
                                                                       131, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 281, 3, 6, 53, 56,
                                                                       137, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 291, 3, 6, 56, 59,
                                                                       143, 149, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 301, 3, 6, 59, 62,
                                                                       149, 155, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 311, 3, 6, 68, 71,
                                                                       161, 167, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 321, 3, 6, 71, 74,
                                                                       167, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 331, 3, 6, 74, 77,
                                                                       173, 179, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 341, 3, 6, 77, 80,
                                                                       179, 185, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 351, 3, 6, 80, 83,
                                                                       185, 191, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 361, 3, 6, 83, 86,
                                                                       191, 197, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 371, 3, 6, 86, 89,
                                                                       197, 203, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 381, 3, 6, 89, 92,
                                                                       203, 209, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 391, 3, 6, 92, 95,
                                                                       209, 215, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 401, 3, 6, 101,
                                                                       107, 221, 231, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 416, 3, 6, 107,
                                                                       113, 231, 241, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 431, 3, 6, 113,
                                                                       119, 241, 251, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 446, 3, 6, 119,
                                                                       125, 251, 261, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 461, 3, 6, 125,
                                                                       131, 261, 271, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 476, 3, 6, 131,
                                                                       137, 271, 281, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 491, 3, 6, 137,
                                                                       143, 281, 291, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 506, 3, 6, 143,
                                                                       149, 291, 301, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 521, 3, 6, 161,
                                                                       167, 311, 321, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 536, 3, 6, 167,
                                                                       173, 321, 331, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 551, 3, 6, 173,
                                                                       179, 331, 341, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 566, 3, 6, 179,
                                                                       185, 341, 351, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 581, 3, 6, 185,
                                                                       191, 351, 361, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 596, 3, 6, 191,
                                                                       197, 361, 371, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 611, 3, 6, 197,
                                                                       203, 371, 381, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 626, 3, 6, 203,
                                                                       209, 381, 391, ncols,
                                                                       gamma, p, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 641, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 644, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 647, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 650, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 653, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 656, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 659, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 662, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 665, 0, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 668, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 671, 0, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 674, 0, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 677, 0, 6, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 680, 0, 6, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 683, 0, 6, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 686, 0, 6, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 689, 0, 6, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 692, 0, 6, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 695, 0, 6, 10, 11,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 704, 0, 6, 11, 12,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 713, 0, 6, 12, 13,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 722, 0, 6, 13, 14,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 731, 0, 6, 14, 15,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 740, 0, 6, 15, 16,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 749, 0, 6, 16, 17,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 758, 0, 6, 17, 18,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 767, 0, 6, 18, 19,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 776, 0, 6, 19, 20,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 785, 0, 6, 23, 24,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 794, 0, 6, 24, 25,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 803, 0, 6, 25, 26,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 812, 0, 6, 26, 27,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 821, 0, 6, 27, 28,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 830, 0, 6, 28, 29,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 839, 0, 6, 29, 30,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 848, 0, 6, 30, 31,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 857, 0, 6, 31, 32,
                                                                       92, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 866, 0, 6, 32, 33,
                                                                       95, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 875, 0, 3, 6, 35,
                                                                       38, 101, 107, 695, 704,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 893, 0, 3, 6, 38,
                                                                       41, 107, 113, 704, 713,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 911, 0, 3, 6, 41,
                                                                       44, 113, 119, 713, 722,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 929, 0, 3, 6, 44,
                                                                       47, 119, 125, 722, 731,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 947, 0, 3, 6, 47,
                                                                       50, 125, 131, 731, 740,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 965, 0, 3, 6, 50,
                                                                       53, 131, 137, 740, 749,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 983, 0, 3, 6, 53,
                                                                       56, 137, 143, 749, 758,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 1001, 0, 3, 6, 56,
                                                                       59, 143, 149, 758, 767,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 1019, 0, 3, 6, 59,
                                                                       62, 149, 155, 767, 776,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 1037, 0, 3, 6, 68,
                                                                       71, 161, 167, 785, 794,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 1055, 0, 3, 6, 71,
                                                                       74, 167, 173, 794, 803,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 1073, 0, 3, 6, 74,
                                                                       77, 173, 179, 803, 812,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 1091, 0, 3, 6, 77,
                                                                       80, 179, 185, 812, 821,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 1109, 0, 3, 6, 80,
                                                                       83, 185, 191, 821, 830,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 1127, 0, 3, 6, 83,
                                                                       86, 191, 197, 830, 839,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 1145, 0, 3, 6, 86,
                                                                       89, 197, 203, 839, 848,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 1163, 0, 3, 6, 89,
                                                                       92, 203, 209, 848, 857,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 1181, 0, 3, 6, 92,
                                                                       95, 209, 215, 857, 866,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1199, 0, 3, 6,
                                                                       101, 107, 221, 231, 875,
                                                                       893, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1229, 0, 3, 6,
                                                                       107, 113, 231, 241, 893,
                                                                       911, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1259, 0, 3, 6,
                                                                       113, 119, 241, 251, 911,
                                                                       929, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1289, 0, 3, 6,
                                                                       119, 125, 251, 261, 929,
                                                                       947, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1319, 0, 3, 6,
                                                                       125, 131, 261, 271, 947,
                                                                       965, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1349, 0, 3, 6,
                                                                       131, 137, 271, 281, 965,
                                                                       983, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1379, 0, 3, 6,
                                                                       137, 143, 281, 291, 983,
                                                                       1001, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1409, 0, 3, 6,
                                                                       143, 149, 291, 301, 1001,
                                                                       1019, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1439, 0, 3, 6,
                                                                       161, 167, 311, 321, 1037,
                                                                       1055, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1469, 0, 3, 6,
                                                                       167, 173, 321, 331, 1055,
                                                                       1073, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1499, 0, 3, 6,
                                                                       173, 179, 331, 341, 1073,
                                                                       1091, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1529, 0, 3, 6,
                                                                       179, 185, 341, 351, 1091,
                                                                       1109, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1559, 0, 3, 6,
                                                                       185, 191, 351, 361, 1109,
                                                                       1127, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1589, 0, 3, 6,
                                                                       191, 197, 361, 371, 1127,
                                                                       1145, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1619, 0, 3, 6,
                                                                       197, 203, 371, 381, 1145,
                                                                       1163, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1649, 0, 3, 6,
                                                                       203, 209, 381, 391, 1163,
                                                                       1181, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1679, 0, 3, 6,
                                                                       221, 231, 401, 416, 1199,
                                                                       1229, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1724, 0, 3, 6,
                                                                       231, 241, 416, 431, 1229,
                                                                       1259, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1769, 0, 3, 6,
                                                                       241, 251, 431, 446, 1259,
                                                                       1289, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1814, 0, 3, 6,
                                                                       251, 261, 446, 461, 1289,
                                                                       1319, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1859, 0, 3, 6,
                                                                       261, 271, 461, 476, 1319,
                                                                       1349, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1904, 0, 3, 6,
                                                                       271, 281, 476, 491, 1349,
                                                                       1379, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1949, 0, 3, 6,
                                                                       281, 291, 491, 506, 1379,
                                                                       1409, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1994, 0, 3, 6,
                                                                       311, 321, 521, 536, 1439,
                                                                       1469, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 2039, 0, 3, 6,
                                                                       321, 331, 536, 551, 1469,
                                                                       1499, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 2084, 0, 3, 6,
                                                                       331, 341, 551, 566, 1499,
                                                                       1529, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 2129, 0, 3, 6,
                                                                       341, 351, 566, 581, 1529,
                                                                       1559, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 2174, 0, 3, 6,
                                                                       351, 361, 581, 596, 1559,
                                                                       1589, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 2219, 0, 3, 6,
                                                                       361, 371, 596, 611, 1589,
                                                                       1619, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 2264, 0, 3, 6,
                                                                       371, 381, 611, 626, 1619,
                                                                       1649, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2309, 6, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2312, 6, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2315, 6, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2318, 6, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2321, 6, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2324, 6, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2327, 6, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2330, 6, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2333, 6, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2336, 6, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2339, 6, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2342, 6, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2345, 6, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2348, 6, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2351, 6, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2354, 6, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2357, 6, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2360, 6, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2363, 6, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2366, 6, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2369, 6, 12, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2378, 6, 13, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2387, 6, 14, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2396, 6, 15, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2405, 6, 16, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2414, 6, 17, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2423, 6, 18, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2432, 6, 19, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2441, 6, 20, 65,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2450, 6, 25, 74,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2459, 6, 26, 77,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2468, 6, 27, 80,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2477, 6, 28, 83,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2486, 6, 29, 86,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2495, 6, 30, 89,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2504, 6, 31, 92,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2513, 6, 32, 95,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2522, 6, 33, 98,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2531, 6, 41, 113,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2549, 6, 44, 119,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2567, 6, 47, 125,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2585, 6, 50, 131,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2603, 6, 53, 137,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2621, 6, 56, 143,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2639, 6, 59, 149,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2657, 6, 62, 155,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2675, 6, 74, 173,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2693, 6, 77, 179,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2711, 6, 80, 185,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2729, 6, 83, 191,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2747, 6, 86, 197,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2765, 6, 89, 203,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2783, 6, 92, 209,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2801, 6, 95, 215,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2819, 6, 113, 241,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2849, 6, 119, 251,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2879, 6, 125, 261,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2909, 6, 131, 271,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2939, 6, 137, 281,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2969, 6, 143, 291,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2999, 6, 149, 301,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3029, 6, 173, 331,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3059, 6, 179, 341,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3089, 6, 185, 351,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3119, 6, 191, 361,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3149, 6, 197, 371,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3179, 6, 203, 381,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3209, 6, 209, 391,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3239, 6, 241, 431,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3284, 6, 251, 446,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3329, 6, 261, 461,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3374, 6, 271, 476,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3419, 6, 281, 491,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3464, 6, 291, 506,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3509, 6, 331, 551,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3554, 6, 341, 566,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3599, 6, 351, 581,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3644, 6, 361, 596,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3689, 6, 371, 611,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3734, 6, 381, 626,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3779, 6, 12, 641,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3788, 6, 13, 644,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3797, 6, 14, 647,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3806, 6, 15, 650,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3815, 6, 16, 653,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3824, 6, 17, 656,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3833, 6, 18, 659,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3842, 6, 19, 662,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3851, 6, 20, 665,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3860, 6, 25, 668,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3869, 6, 26, 671,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3878, 6, 27, 674,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3887, 6, 28, 677,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3896, 6, 29, 680,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3905, 6, 30, 683,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3914, 6, 31, 686,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3923, 6, 32, 689,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3932, 6, 33, 692,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3941, 6, 41, 641,
                                                                       713, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3968, 6, 44, 644,
                                                                       722, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3995, 6, 47, 647,
                                                                       731, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 4022, 6, 50, 650,
                                                                       740, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 4049, 6, 53, 653,
                                                                       749, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 4076, 6, 56, 656,
                                                                       758, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 4103, 6, 59, 659,
                                                                       767, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 4130, 6, 62, 662,
                                                                       776, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 4157, 6, 74, 668,
                                                                       803, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 4184, 6, 77, 671,
                                                                       812, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 4211, 6, 80, 674,
                                                                       821, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 4238, 6, 83, 677,
                                                                       830, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 4265, 6, 86, 680,
                                                                       839, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 4292, 6, 89, 683,
                                                                       848, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 4319, 6, 92, 686,
                                                                       857, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 4346, 6, 95, 689,
                                                                       866, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4373, 3, 6, 113,
                                                                       3941, 713, 3968, 911,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4427, 3, 6, 119,
                                                                       3968, 722, 3995, 929,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4481, 3, 6, 125,
                                                                       3995, 731, 4022, 947,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4535, 3, 6, 131,
                                                                       4022, 740, 4049, 965,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4589, 3, 6, 137,
                                                                       4049, 749, 4076, 983,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4643, 3, 6, 143,
                                                                       4076, 758, 4103, 1001,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4697, 3, 6, 149,
                                                                       4103, 767, 4130, 1019,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4751, 3, 6, 173,
                                                                       4157, 803, 4184, 1073,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4805, 3, 6, 179,
                                                                       4184, 812, 4211, 1091,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4859, 3, 6, 185,
                                                                       4211, 821, 4238, 1109,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4913, 3, 6, 191,
                                                                       4238, 830, 4265, 1127,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4967, 3, 6, 197,
                                                                       4265, 839, 4292, 1145,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 5021, 3, 6, 203,
                                                                       4292, 848, 4319, 1163,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 5075, 3, 6, 209,
                                                                       4319, 857, 4346, 1181,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5129, 3, 6, 241,
                                                                       4373, 911, 4427, 1259,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5219, 3, 6, 251,
                                                                       4427, 929, 4481, 1289,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5309, 3, 6, 261,
                                                                       4481, 947, 4535, 1319,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5399, 3, 6, 271,
                                                                       4535, 965, 4589, 1349,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5489, 3, 6, 281,
                                                                       4589, 983, 4643, 1379,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5579, 3, 6, 291,
                                                                       4643, 1001, 4697, 1409,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5669, 3, 6, 331,
                                                                       4751, 1073, 4805, 1499,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5759, 3, 6, 341,
                                                                       4805, 1091, 4859, 1529,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5849, 3, 6, 351,
                                                                       4859, 1109, 4913, 1559,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5939, 3, 6, 361,
                                                                       4913, 1127, 4967, 1589,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 6029, 3, 6, 371,
                                                                       4967, 1145, 5021, 1619,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 6119, 3, 6, 381,
                                                                       5021, 1163, 5075, 1649,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 6209, 3, 6, 431,
                                                                       5129, 1259, 5219, 1769,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 6344, 3, 6, 446,
                                                                       5219, 1289, 5309, 1814,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 6479, 3, 6, 461,
                                                                       5309, 1319, 5399, 1859,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 6614, 3, 6, 476,
                                                                       5399, 1349, 5489, 1904,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 6749, 3, 6, 491,
                                                                       5489, 1379, 5579, 1949,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 6884, 3, 6, 551,
                                                                       5669, 1499, 5759, 2084,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 7019, 3, 6, 566,
                                                                       5759, 1529, 5849, 2129,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 7154, 3, 6, 581,
                                                                       5849, 1559, 5939, 2174,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 7289, 3, 6, 596,
                                                                       5939, 1589, 6029, 2219,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 7424, 3, 6, 611,
                                                                       6029, 1619, 6119, 2264,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7559, 6, 10, 11,
                                                                       2309, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7565, 6, 11, 12,
                                                                       2312, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7571, 6, 12, 13,
                                                                       2315, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7577, 6, 13, 14,
                                                                       2318, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7583, 6, 14, 15,
                                                                       2321, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7589, 6, 15, 16,
                                                                       2324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7595, 6, 16, 17,
                                                                       2327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7601, 6, 17, 18,
                                                                       2330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7607, 6, 18, 19,
                                                                       2333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7613, 6, 19, 20,
                                                                       2336, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7619, 6, 23, 24,
                                                                       2339, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7625, 6, 24, 25,
                                                                       2342, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7631, 6, 25, 26,
                                                                       2345, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7637, 6, 26, 27,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7643, 6, 27, 28,
                                                                       2351, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7649, 6, 28, 29,
                                                                       2354, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7655, 6, 29, 30,
                                                                       2357, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7661, 6, 30, 31,
                                                                       2360, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7667, 6, 31, 32,
                                                                       2363, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7673, 6, 32, 33,
                                                                       2366, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7679, 3, 6, 7559,
                                                                       2309, 7565, 2369, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7697, 3, 6, 7565,
                                                                       2312, 7571, 2378, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7715, 3, 6, 7571,
                                                                       2315, 7577, 2387, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7733, 3, 6, 7577,
                                                                       2318, 7583, 2396, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7751, 3, 6, 7583,
                                                                       2321, 7589, 2405, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7769, 3, 6, 7589,
                                                                       2324, 7595, 2414, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7787, 3, 6, 7595,
                                                                       2327, 7601, 2423, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7805, 3, 6, 7601,
                                                                       2330, 7607, 2432, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7823, 3, 6, 7607,
                                                                       2333, 7613, 2441, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7841, 3, 6, 7619,
                                                                       2339, 7625, 2450, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7859, 3, 6, 7625,
                                                                       2342, 7631, 2459, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7877, 3, 6, 7631,
                                                                       2345, 7637, 2468, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7895, 3, 6, 7637,
                                                                       2348, 7643, 2477, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7913, 3, 6, 7643,
                                                                       2351, 7649, 2486, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7931, 3, 6, 7649,
                                                                       2354, 7655, 2495, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7949, 3, 6, 7655,
                                                                       2357, 7661, 2504, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7967, 3, 6, 7661,
                                                                       2360, 7667, 2513, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7985, 3, 6, 7667,
                                                                       2363, 7673, 2522, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8003, 3, 6, 7679,
                                                                       2369, 7697, 101, 107,
                                                                       2531, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8039, 3, 6, 7697,
                                                                       2378, 7715, 107, 113,
                                                                       2549, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8075, 3, 6, 7715,
                                                                       2387, 7733, 113, 119,
                                                                       2567, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8111, 3, 6, 7733,
                                                                       2396, 7751, 119, 125,
                                                                       2585, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8147, 3, 6, 7751,
                                                                       2405, 7769, 125, 131,
                                                                       2603, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8183, 3, 6, 7769,
                                                                       2414, 7787, 131, 137,
                                                                       2621, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8219, 3, 6, 7787,
                                                                       2423, 7805, 137, 143,
                                                                       2639, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8255, 3, 6, 7805,
                                                                       2432, 7823, 143, 149,
                                                                       2657, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8291, 3, 6, 7841,
                                                                       2450, 7859, 161, 167,
                                                                       2675, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8327, 3, 6, 7859,
                                                                       2459, 7877, 167, 173,
                                                                       2693, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8363, 3, 6, 7877,
                                                                       2468, 7895, 173, 179,
                                                                       2711, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8399, 3, 6, 7895,
                                                                       2477, 7913, 179, 185,
                                                                       2729, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8435, 3, 6, 7913,
                                                                       2486, 7931, 185, 191,
                                                                       2747, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8471, 3, 6, 7931,
                                                                       2495, 7949, 191, 197,
                                                                       2765, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8507, 3, 6, 7949,
                                                                       2504, 7967, 197, 203,
                                                                       2783, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8543, 3, 6, 7967,
                                                                       2513, 7985, 203, 209,
                                                                       2801, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8579, 3, 6, 8003,
                                                                       2531, 8039, 221, 231,
                                                                       2819, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8639, 3, 6, 8039,
                                                                       2549, 8075, 231, 241,
                                                                       2849, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8699, 3, 6, 8075,
                                                                       2567, 8111, 241, 251,
                                                                       2879, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8759, 3, 6, 8111,
                                                                       2585, 8147, 251, 261,
                                                                       2909, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8819, 3, 6, 8147,
                                                                       2603, 8183, 261, 271,
                                                                       2939, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8879, 3, 6, 8183,
                                                                       2621, 8219, 271, 281,
                                                                       2969, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8939, 3, 6, 8219,
                                                                       2639, 8255, 281, 291,
                                                                       2999, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8999, 3, 6, 8291,
                                                                       2675, 8327, 311, 321,
                                                                       3029, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9059, 3, 6, 8327,
                                                                       2693, 8363, 321, 331,
                                                                       3059, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9119, 3, 6, 8363,
                                                                       2711, 8399, 331, 341,
                                                                       3089, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9179, 3, 6, 8399,
                                                                       2729, 8435, 341, 351,
                                                                       3119, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9239, 3, 6, 8435,
                                                                       2747, 8471, 351, 361,
                                                                       3149, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9299, 3, 6, 8471,
                                                                       2765, 8507, 361, 371,
                                                                       3179, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9359, 3, 6, 8507,
                                                                       2783, 8543, 371, 381,
                                                                       3209, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9419, 3, 6, 8579,
                                                                       2819, 8639, 401, 416,
                                                                       3239, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9509, 3, 6, 8639,
                                                                       2849, 8699, 416, 431,
                                                                       3284, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9599, 3, 6, 8699,
                                                                       2879, 8759, 431, 446,
                                                                       3329, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9689, 3, 6, 8759,
                                                                       2909, 8819, 446, 461,
                                                                       3374, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9779, 3, 6, 8819,
                                                                       2939, 8879, 461, 476,
                                                                       3419, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9869, 3, 6, 8879,
                                                                       2969, 8939, 476, 491,
                                                                       3464, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9959, 3, 6, 8999,
                                                                       3029, 9059, 521, 536,
                                                                       3509, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10049, 3, 6, 9059,
                                                                       3059, 9119, 536, 551,
                                                                       3554, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10139, 3, 6, 9119,
                                                                       3089, 9179, 551, 566,
                                                                       3599, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10229, 3, 6, 9179,
                                                                       3119, 9239, 566, 581,
                                                                       3644, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10319, 3, 6, 9239,
                                                                       3149, 9299, 581, 596,
                                                                       3689, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10409, 3, 6, 9299,
                                                                       3179, 9359, 596, 611,
                                                                       3734, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10499, 0, 6, 7559,
                                                                       2309, 7565, 3779, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10517, 0, 6, 7565,
                                                                       2312, 7571, 3788, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10535, 0, 6, 7571,
                                                                       2315, 7577, 3797, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10553, 0, 6, 7577,
                                                                       2318, 7583, 3806, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10571, 0, 6, 7583,
                                                                       2321, 7589, 3815, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10589, 0, 6, 7589,
                                                                       2324, 7595, 3824, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10607, 0, 6, 7595,
                                                                       2327, 7601, 3833, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10625, 0, 6, 7601,
                                                                       2330, 7607, 3842, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10643, 0, 6, 7607,
                                                                       2333, 7613, 3851, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10661, 0, 6, 7619,
                                                                       2339, 7625, 3860, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10679, 0, 6, 7625,
                                                                       2342, 7631, 3869, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10697, 0, 6, 7631,
                                                                       2345, 7637, 3878, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10715, 0, 6, 7637,
                                                                       2348, 7643, 3887, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10733, 0, 6, 7643,
                                                                       2351, 7649, 3896, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10751, 0, 6, 7649,
                                                                       2354, 7655, 3905, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10769, 0, 6, 7655,
                                                                       2357, 7661, 3914, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10787, 0, 6, 7661,
                                                                       2360, 7667, 3923, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 10805, 0, 6, 7667,
                                                                       2363, 7673, 3932, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 10823, 0, 3, 6,
                                                                       7679, 2369, 7697, 10499,
                                                                       3779, 10517, 695, 704,
                                                                       3941, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 10877, 0, 3, 6,
                                                                       7697, 2378, 7715, 10517,
                                                                       3788, 10535, 704, 713,
                                                                       3968, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 10931, 0, 3, 6,
                                                                       7715, 2387, 7733, 10535,
                                                                       3797, 10553, 713, 722,
                                                                       3995, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 10985, 0, 3, 6,
                                                                       7733, 2396, 7751, 10553,
                                                                       3806, 10571, 722, 731,
                                                                       4022, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 11039, 0, 3, 6,
                                                                       7751, 2405, 7769, 10571,
                                                                       3815, 10589, 731, 740,
                                                                       4049, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 11093, 0, 3, 6,
                                                                       7769, 2414, 7787, 10589,
                                                                       3824, 10607, 740, 749,
                                                                       4076, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 11147, 0, 3, 6,
                                                                       7787, 2423, 7805, 10607,
                                                                       3833, 10625, 749, 758,
                                                                       4103, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 11201, 0, 3, 6,
                                                                       7805, 2432, 7823, 10625,
                                                                       3842, 10643, 758, 767,
                                                                       4130, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 11255, 0, 3, 6,
                                                                       7841, 2450, 7859, 10661,
                                                                       3860, 10679, 785, 794,
                                                                       4157, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 11309, 0, 3, 6,
                                                                       7859, 2459, 7877, 10679,
                                                                       3869, 10697, 794, 803,
                                                                       4184, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 11363, 0, 3, 6,
                                                                       7877, 2468, 7895, 10697,
                                                                       3878, 10715, 803, 812,
                                                                       4211, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 11417, 0, 3, 6,
                                                                       7895, 2477, 7913, 10715,
                                                                       3887, 10733, 812, 821,
                                                                       4238, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 11471, 0, 3, 6,
                                                                       7913, 2486, 7931, 10733,
                                                                       3896, 10751, 821, 830,
                                                                       4265, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 11525, 0, 3, 6,
                                                                       7931, 2495, 7949, 10751,
                                                                       3905, 10769, 830, 839,
                                                                       4292, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 11579, 0, 3, 6,
                                                                       7949, 2504, 7967, 10769,
                                                                       3914, 10787, 839, 848,
                                                                       4319, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 11633, 0, 3, 6,
                                                                       7967, 2513, 7985, 10787,
                                                                       3923, 10805, 848, 857,
                                                                       4346, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 11687, 0, 3, 6,
                                                                       8003, 2531, 8039, 10823,
                                                                       3941, 10877, 875, 893,
                                                                       4373, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 11795, 0, 3, 6,
                                                                       8039, 2549, 8075, 10877,
                                                                       3968, 10931, 893, 911,
                                                                       4427, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 11903, 0, 3, 6,
                                                                       8075, 2567, 8111, 10931,
                                                                       3995, 10985, 911, 929,
                                                                       4481, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 12011, 0, 3, 6,
                                                                       8111, 2585, 8147, 10985,
                                                                       4022, 11039, 929, 947,
                                                                       4535, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 12119, 0, 3, 6,
                                                                       8147, 2603, 8183, 11039,
                                                                       4049, 11093, 947, 965,
                                                                       4589, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 12227, 0, 3, 6,
                                                                       8183, 2621, 8219, 11093,
                                                                       4076, 11147, 965, 983,
                                                                       4643, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 12335, 0, 3, 6,
                                                                       8219, 2639, 8255, 11147,
                                                                       4103, 11201, 983, 1001,
                                                                       4697, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 12443, 0, 3, 6,
                                                                       8291, 2675, 8327, 11255,
                                                                       4157, 11309, 1037, 1055,
                                                                       4751, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 12551, 0, 3, 6,
                                                                       8327, 2693, 8363, 11309,
                                                                       4184, 11363, 1055, 1073,
                                                                       4805, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 12659, 0, 3, 6,
                                                                       8363, 2711, 8399, 11363,
                                                                       4211, 11417, 1073, 1091,
                                                                       4859, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 12767, 0, 3, 6,
                                                                       8399, 2729, 8435, 11417,
                                                                       4238, 11471, 1091, 1109,
                                                                       4913, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 12875, 0, 3, 6,
                                                                       8435, 2747, 8471, 11471,
                                                                       4265, 11525, 1109, 1127,
                                                                       4967, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 12983, 0, 3, 6,
                                                                       8471, 2765, 8507, 11525,
                                                                       4292, 11579, 1127, 1145,
                                                                       5021, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 13091, 0, 3, 6,
                                                                       8507, 2783, 8543, 11579,
                                                                       4319, 11633, 1145, 1163,
                                                                       5075, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 13199, 0, 3, 6,
                                                                       8579, 2819, 8639, 10823,
                                                                       10877, 11687, 4373, 11795,
                                                                       1199, 1229, 5129, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 13379, 0, 3, 6,
                                                                       8639, 2849, 8699, 10877,
                                                                       10931, 11795, 4427, 11903,
                                                                       1229, 1259, 5219, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 13559, 0, 3, 6,
                                                                       8699, 2879, 8759, 10931,
                                                                       10985, 11903, 4481, 12011,
                                                                       1259, 1289, 5309, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 13739, 0, 3, 6,
                                                                       8759, 2909, 8819, 10985,
                                                                       11039, 12011, 4535, 12119,
                                                                       1289, 1319, 5399, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 13919, 0, 3, 6,
                                                                       8819, 2939, 8879, 11039,
                                                                       11093, 12119, 4589, 12227,
                                                                       1319, 1349, 5489, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 14099, 0, 3, 6,
                                                                       8879, 2969, 8939, 11093,
                                                                       11147, 12227, 4643, 12335,
                                                                       1349, 1379, 5579, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 14279, 0, 3, 6,
                                                                       8999, 3029, 9059, 11255,
                                                                       11309, 12443, 4751, 12551,
                                                                       1439, 1469, 5669, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 14459, 0, 3, 6,
                                                                       9059, 3059, 9119, 11309,
                                                                       11363, 12551, 4805, 12659,
                                                                       1469, 1499, 5759, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 14639, 0, 3, 6,
                                                                       9119, 3089, 9179, 11363,
                                                                       11417, 12659, 4859, 12767,
                                                                       1499, 1529, 5849, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 14819, 0, 3, 6,
                                                                       9179, 3119, 9239, 11417,
                                                                       11471, 12767, 4913, 12875,
                                                                       1529, 1559, 5939, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 14999, 0, 3, 6,
                                                                       9239, 3149, 9299, 11471,
                                                                       11525, 12875, 4967, 12983,
                                                                       1559, 1589, 6029, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 15179, 0, 3, 6,
                                                                       9299, 3179, 9359, 11525,
                                                                       11579, 12983, 5021, 13091,
                                                                       1589, 1619, 6119, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 15359, 0, 3, 6,
                                                                       9419, 3239, 9509, 11687,
                                                                       11795, 13199, 5129, 13379,
                                                                       1679, 1724, 6209, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 15629, 0, 3, 6,
                                                                       9509, 3284, 9599, 11795,
                                                                       11903, 13379, 5219, 13559,
                                                                       1724, 1769, 6344, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 15899, 0, 3, 6,
                                                                       9599, 3329, 9689, 11903,
                                                                       12011, 13559, 5309, 13739,
                                                                       1769, 1814, 6479, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 16169, 0, 3, 6,
                                                                       9689, 3374, 9779, 12011,
                                                                       12119, 13739, 5399, 13919,
                                                                       1814, 1859, 6614, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 16439, 0, 3, 6,
                                                                       9779, 3419, 9869, 12119,
                                                                       12227, 13919, 5489, 14099,
                                                                       1859, 1904, 6749, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 16709, 0, 3, 6,
                                                                       9959, 3509, 10049, 12443,
                                                                       12551, 14279, 5669, 14459,
                                                                       1994, 2039, 6884, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 16979, 0, 3, 6,
                                                                       10049, 3554, 10139, 12551,
                                                                       12659, 14459, 5759, 14639,
                                                                       2039, 2084, 7019, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 17249, 0, 3, 6,
                                                                       10139, 3599, 10229, 12659,
                                                                       12767, 14639, 5849, 14819,
                                                                       2084, 2129, 7154, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 17519, 0, 3, 6,
                                                                       10229, 3644, 10319, 12767,
                                                                       12875, 14819, 5939, 14999,
                                                                       2129, 2174, 7289, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 17789, 0, 3, 6,
                                                                       10319, 3689, 10409, 12875,
                                                                       12983, 14999, 6029, 15179,
                                                                       2174, 2219, 7424, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18059, 6, 2309,
                                                                       2312, 7571, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18069, 6, 2312,
                                                                       2315, 7577, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18079, 6, 2315,
                                                                       2318, 7583, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18089, 6, 2318,
                                                                       2321, 7589, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18099, 6, 2321,
                                                                       2324, 7595, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18109, 6, 2324,
                                                                       2327, 7601, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18119, 6, 2327,
                                                                       2330, 7607, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18129, 6, 2330,
                                                                       2333, 7613, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18139, 6, 2339,
                                                                       2342, 7631, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18149, 6, 2342,
                                                                       2345, 7637, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18159, 6, 2345,
                                                                       2348, 7643, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18169, 6, 2348,
                                                                       2351, 7649, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18179, 6, 2351,
                                                                       2354, 7655, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18189, 6, 2354,
                                                                       2357, 7661, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18199, 6, 2357,
                                                                       2360, 7667, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18209, 6, 2360,
                                                                       2363, 7673, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18219, 3, 6,
                                                                       18059, 7571, 18069, 7715,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18249, 3, 6,
                                                                       18069, 7577, 18079, 7733,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18279, 3, 6,
                                                                       18079, 7583, 18089, 7751,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18309, 3, 6,
                                                                       18089, 7589, 18099, 7769,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18339, 3, 6,
                                                                       18099, 7595, 18109, 7787,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18369, 3, 6,
                                                                       18109, 7601, 18119, 7805,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18399, 3, 6,
                                                                       18119, 7607, 18129, 7823,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18429, 3, 6,
                                                                       18139, 7631, 18149, 7877,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18459, 3, 6,
                                                                       18149, 7637, 18159, 7895,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18489, 3, 6,
                                                                       18159, 7643, 18169, 7913,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18519, 3, 6,
                                                                       18169, 7649, 18179, 7931,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18549, 3, 6,
                                                                       18179, 7655, 18189, 7949,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18579, 3, 6,
                                                                       18189, 7661, 18199, 7967,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18609, 3, 6,
                                                                       18199, 7667, 18209, 7985,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18639, 3, 6,
                                                                       18219, 7715, 18249, 2531,
                                                                       2549, 8075, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18699, 3, 6,
                                                                       18249, 7733, 18279, 2549,
                                                                       2567, 8111, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18759, 3, 6,
                                                                       18279, 7751, 18309, 2567,
                                                                       2585, 8147, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18819, 3, 6,
                                                                       18309, 7769, 18339, 2585,
                                                                       2603, 8183, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18879, 3, 6,
                                                                       18339, 7787, 18369, 2603,
                                                                       2621, 8219, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18939, 3, 6,
                                                                       18369, 7805, 18399, 2621,
                                                                       2639, 8255, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18999, 3, 6,
                                                                       18429, 7877, 18459, 2675,
                                                                       2693, 8363, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19059, 3, 6,
                                                                       18459, 7895, 18489, 2693,
                                                                       2711, 8399, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19119, 3, 6,
                                                                       18489, 7913, 18519, 2711,
                                                                       2729, 8435, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19179, 3, 6,
                                                                       18519, 7931, 18549, 2729,
                                                                       2747, 8471, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19239, 3, 6,
                                                                       18549, 7949, 18579, 2747,
                                                                       2765, 8507, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19299, 3, 6,
                                                                       18579, 7967, 18609, 2765,
                                                                       2783, 8543, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19359, 3, 6,
                                                                       18639, 8075, 18699, 2819,
                                                                       2849, 8699, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19459, 3, 6,
                                                                       18699, 8111, 18759, 2849,
                                                                       2879, 8759, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19559, 3, 6,
                                                                       18759, 8147, 18819, 2879,
                                                                       2909, 8819, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19659, 3, 6,
                                                                       18819, 8183, 18879, 2909,
                                                                       2939, 8879, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19759, 3, 6,
                                                                       18879, 8219, 18939, 2939,
                                                                       2969, 8939, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19859, 3, 6,
                                                                       18999, 8363, 19059, 3029,
                                                                       3059, 9119, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19959, 3, 6,
                                                                       19059, 8399, 19119, 3059,
                                                                       3089, 9179, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20059, 3, 6,
                                                                       19119, 8435, 19179, 3089,
                                                                       3119, 9239, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20159, 3, 6,
                                                                       19179, 8471, 19239, 3119,
                                                                       3149, 9299, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20259, 3, 6,
                                                                       19239, 8507, 19299, 3149,
                                                                       3179, 9359, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20359, 3, 6,
                                                                       19359, 8699, 19459, 3239,
                                                                       3284, 9599, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20509, 3, 6,
                                                                       19459, 8759, 19559, 3284,
                                                                       3329, 9689, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20659, 3, 6,
                                                                       19559, 8819, 19659, 3329,
                                                                       3374, 9779, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20809, 3, 6,
                                                                       19659, 8879, 19759, 3374,
                                                                       3419, 9869, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20959, 3, 6,
                                                                       19859, 9119, 19959, 3509,
                                                                       3554, 10139, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21109, 3, 6,
                                                                       19959, 9179, 20059, 3554,
                                                                       3599, 10229, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21259, 3, 6,
                                                                       20059, 9239, 20159, 3599,
                                                                       3644, 10319, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21409, 3, 6,
                                                                       20159, 9299, 20259, 3644,
                                                                       3689, 10409, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21559, 0, 6,
                                                                       18059, 7571, 18069, 3779,
                                                                       3788, 10535, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21589, 0, 6,
                                                                       18069, 7577, 18079, 3788,
                                                                       3797, 10553, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21619, 0, 6,
                                                                       18079, 7583, 18089, 3797,
                                                                       3806, 10571, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21649, 0, 6,
                                                                       18089, 7589, 18099, 3806,
                                                                       3815, 10589, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21679, 0, 6,
                                                                       18099, 7595, 18109, 3815,
                                                                       3824, 10607, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21709, 0, 6,
                                                                       18109, 7601, 18119, 3824,
                                                                       3833, 10625, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21739, 0, 6,
                                                                       18119, 7607, 18129, 3833,
                                                                       3842, 10643, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21769, 0, 6,
                                                                       18139, 7631, 18149, 3860,
                                                                       3869, 10697, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21799, 0, 6,
                                                                       18149, 7637, 18159, 3869,
                                                                       3878, 10715, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21829, 0, 6,
                                                                       18159, 7643, 18169, 3878,
                                                                       3887, 10733, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21859, 0, 6,
                                                                       18169, 7649, 18179, 3887,
                                                                       3896, 10751, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21889, 0, 6,
                                                                       18179, 7655, 18189, 3896,
                                                                       3905, 10769, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21919, 0, 6,
                                                                       18189, 7661, 18199, 3905,
                                                                       3914, 10787, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21949, 0, 6,
                                                                       18199, 7667, 18209, 3914,
                                                                       3923, 10805, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 21979, 0, 3, 6,
                                                                       18219, 7715, 18249, 21559,
                                                                       10535, 21589, 3941, 3968,
                                                                       10931, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 22069, 0, 3, 6,
                                                                       18249, 7733, 18279, 21589,
                                                                       10553, 21619, 3968, 3995,
                                                                       10985, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 22159, 0, 3, 6,
                                                                       18279, 7751, 18309, 21619,
                                                                       10571, 21649, 3995, 4022,
                                                                       11039, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 22249, 0, 3, 6,
                                                                       18309, 7769, 18339, 21649,
                                                                       10589, 21679, 4022, 4049,
                                                                       11093, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 22339, 0, 3, 6,
                                                                       18339, 7787, 18369, 21679,
                                                                       10607, 21709, 4049, 4076,
                                                                       11147, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 22429, 0, 3, 6,
                                                                       18369, 7805, 18399, 21709,
                                                                       10625, 21739, 4076, 4103,
                                                                       11201, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 22519, 0, 3, 6,
                                                                       18429, 7877, 18459, 21769,
                                                                       10697, 21799, 4157, 4184,
                                                                       11363, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 22609, 0, 3, 6,
                                                                       18459, 7895, 18489, 21799,
                                                                       10715, 21829, 4184, 4211,
                                                                       11417, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 22699, 0, 3, 6,
                                                                       18489, 7913, 18519, 21829,
                                                                       10733, 21859, 4211, 4238,
                                                                       11471, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 22789, 0, 3, 6,
                                                                       18519, 7931, 18549, 21859,
                                                                       10751, 21889, 4238, 4265,
                                                                       11525, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 22879, 0, 3, 6,
                                                                       18549, 7949, 18579, 21889,
                                                                       10769, 21919, 4265, 4292,
                                                                       11579, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 22969, 0, 3, 6,
                                                                       18579, 7967, 18609, 21919,
                                                                       10787, 21949, 4292, 4319,
                                                                       11633, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 23059, 0, 3, 6,
                                                                       18639, 8075, 18699, 21979,
                                                                       10931, 22069, 4373, 4427,
                                                                       11903, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 23239, 0, 3, 6,
                                                                       18699, 8111, 18759, 22069,
                                                                       10985, 22159, 4427, 4481,
                                                                       12011, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 23419, 0, 3, 6,
                                                                       18759, 8147, 18819, 22159,
                                                                       11039, 22249, 4481, 4535,
                                                                       12119, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 23599, 0, 3, 6,
                                                                       18819, 8183, 18879, 22249,
                                                                       11093, 22339, 4535, 4589,
                                                                       12227, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 23779, 0, 3, 6,
                                                                       18879, 8219, 18939, 22339,
                                                                       11147, 22429, 4589, 4643,
                                                                       12335, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 23959, 0, 3, 6,
                                                                       18999, 8363, 19059, 22519,
                                                                       11363, 22609, 4751, 4805,
                                                                       12659, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 24139, 0, 3, 6,
                                                                       19059, 8399, 19119, 22609,
                                                                       11417, 22699, 4805, 4859,
                                                                       12767, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 24319, 0, 3, 6,
                                                                       19119, 8435, 19179, 22699,
                                                                       11471, 22789, 4859, 4913,
                                                                       12875, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 24499, 0, 3, 6,
                                                                       19179, 8471, 19239, 22789,
                                                                       11525, 22879, 4913, 4967,
                                                                       12983, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 24679, 0, 3, 6,
                                                                       19239, 8507, 19299, 22879,
                                                                       11579, 22969, 4967, 5021,
                                                                       13091, ncols, gamma, p,
                                                                       q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 24859, 0, 3, 6,
                                                                       19359, 8699, 19459, 21979,
                                                                       22069, 23059, 11903,
                                                                       23239, 5129, 5219, 13559,
                                                                       ncols, gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 25159, 0, 3, 6,
                                                                       19459, 8759, 19559, 22069,
                                                                       22159, 23239, 12011,
                                                                       23419, 5219, 5309, 13739,
                                                                       ncols, gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 25459, 0, 3, 6,
                                                                       19559, 8819, 19659, 22159,
                                                                       22249, 23419, 12119,
                                                                       23599, 5309, 5399, 13919,
                                                                       ncols, gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 25759, 0, 3, 6,
                                                                       19659, 8879, 19759, 22249,
                                                                       22339, 23599, 12227,
                                                                       23779, 5399, 5489, 14099,
                                                                       ncols, gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 26059, 0, 3, 6,
                                                                       19859, 9119, 19959, 22519,
                                                                       22609, 23959, 12659,
                                                                       24139, 5669, 5759, 14639,
                                                                       ncols, gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 26359, 0, 3, 6,
                                                                       19959, 9179, 20059, 22609,
                                                                       22699, 24139, 12767,
                                                                       24319, 5759, 5849, 14819,
                                                                       ncols, gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 26659, 0, 3, 6,
                                                                       20059, 9239, 20159, 22699,
                                                                       22789, 24319, 12875,
                                                                       24499, 5849, 5939, 14999,
                                                                       ncols, gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 26959, 0, 3, 6,
                                                                       20159, 9299, 20259, 22789,
                                                                       22879, 24499, 12983,
                                                                       24679, 5939, 6029, 15179,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 27259, 0, 3, 6,
                                                                       20359, 9599, 20509, 23059,
                                                                       23239, 24859, 13559,
                                                                       25159, 6209, 6344, 15899,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 27709, 0, 3, 6,
                                                                       20509, 9689, 20659, 23239,
                                                                       23419, 25159, 13739,
                                                                       25459, 6344, 6479, 16169,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 28159, 0, 3, 6,
                                                                       20659, 9779, 20809, 23419,
                                                                       23599, 25459, 13919,
                                                                       25759, 6479, 6614, 16439,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 28609, 0, 3, 6,
                                                                       20959, 10139, 21109,
                                                                       23959, 24139, 26059,
                                                                       14639, 26359, 6884, 7019,
                                                                       17249, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 29059, 0, 3, 6,
                                                                       21109, 10229, 21259,
                                                                       24139, 24319, 26359,
                                                                       14819, 26659, 7019, 7154,
                                                                       17519, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 29509, 0, 3, 6,
                                                                       21259, 10319, 21409,
                                                                       24319, 24499, 26659,
                                                                       14999, 26959, 7154, 7289,
                                                                       17789, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 29959, 6, 7559,
                                                                       7565, 18059, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 29974, 6, 7565,
                                                                       7571, 18069, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 29989, 6, 7571,
                                                                       7577, 18079, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30004, 6, 7577,
                                                                       7583, 18089, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30019, 6, 7583,
                                                                       7589, 18099, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30034, 6, 7589,
                                                                       7595, 18109, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30049, 6, 7595,
                                                                       7601, 18119, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30064, 6, 7601,
                                                                       7607, 18129, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30079, 6, 7619,
                                                                       7625, 18139, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30094, 6, 7625,
                                                                       7631, 18149, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30109, 6, 7631,
                                                                       7637, 18159, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30124, 6, 7637,
                                                                       7643, 18169, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30139, 6, 7643,
                                                                       7649, 18179, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30154, 6, 7649,
                                                                       7655, 18189, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30169, 6, 7655,
                                                                       7661, 18199, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30184, 6, 7661,
                                                                       7667, 18209, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30199, 3, 6,
                                                                       29959, 18059, 29974, 7679,
                                                                       7697, 18219, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30244, 3, 6,
                                                                       29974, 18069, 29989, 7697,
                                                                       7715, 18249, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30289, 3, 6,
                                                                       29989, 18079, 30004, 7715,
                                                                       7733, 18279, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30334, 3, 6,
                                                                       30004, 18089, 30019, 7733,
                                                                       7751, 18309, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30379, 3, 6,
                                                                       30019, 18099, 30034, 7751,
                                                                       7769, 18339, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30424, 3, 6,
                                                                       30034, 18109, 30049, 7769,
                                                                       7787, 18369, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30469, 3, 6,
                                                                       30049, 18119, 30064, 7787,
                                                                       7805, 18399, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30514, 3, 6,
                                                                       30079, 18139, 30094, 7841,
                                                                       7859, 18429, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30559, 3, 6,
                                                                       30094, 18149, 30109, 7859,
                                                                       7877, 18459, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30604, 3, 6,
                                                                       30109, 18159, 30124, 7877,
                                                                       7895, 18489, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30649, 3, 6,
                                                                       30124, 18169, 30139, 7895,
                                                                       7913, 18519, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30694, 3, 6,
                                                                       30139, 18179, 30154, 7913,
                                                                       7931, 18549, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30739, 3, 6,
                                                                       30154, 18189, 30169, 7931,
                                                                       7949, 18579, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30784, 3, 6,
                                                                       30169, 18199, 30184, 7949,
                                                                       7967, 18609, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 30829, 3, 6,
                                                                       30199, 18219, 30244, 8003,
                                                                       8039, 18639, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 30919, 3, 6,
                                                                       30244, 18249, 30289, 8039,
                                                                       8075, 18699, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 31009, 3, 6,
                                                                       30289, 18279, 30334, 8075,
                                                                       8111, 18759, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 31099, 3, 6,
                                                                       30334, 18309, 30379, 8111,
                                                                       8147, 18819, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 31189, 3, 6,
                                                                       30379, 18339, 30424, 8147,
                                                                       8183, 18879, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 31279, 3, 6,
                                                                       30424, 18369, 30469, 8183,
                                                                       8219, 18939, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 31369, 3, 6,
                                                                       30514, 18429, 30559, 8291,
                                                                       8327, 18999, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 31459, 3, 6,
                                                                       30559, 18459, 30604, 8327,
                                                                       8363, 19059, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 31549, 3, 6,
                                                                       30604, 18489, 30649, 8363,
                                                                       8399, 19119, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 31639, 3, 6,
                                                                       30649, 18519, 30694, 8399,
                                                                       8435, 19179, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 31729, 3, 6,
                                                                       30694, 18549, 30739, 8435,
                                                                       8471, 19239, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 31819, 3, 6,
                                                                       30739, 18579, 30784, 8471,
                                                                       8507, 19299, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 31909, 3, 6,
                                                                       30829, 18639, 30919, 8579,
                                                                       8639, 19359, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 32059, 3, 6,
                                                                       30919, 18699, 31009, 8639,
                                                                       8699, 19459, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 32209, 3, 6,
                                                                       31009, 18759, 31099, 8699,
                                                                       8759, 19559, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 32359, 3, 6,
                                                                       31099, 18819, 31189, 8759,
                                                                       8819, 19659, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 32509, 3, 6,
                                                                       31189, 18879, 31279, 8819,
                                                                       8879, 19759, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 32659, 3, 6,
                                                                       31369, 18999, 31459, 8999,
                                                                       9059, 19859, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 32809, 3, 6,
                                                                       31459, 19059, 31549, 9059,
                                                                       9119, 19959, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 32959, 3, 6,
                                                                       31549, 19119, 31639, 9119,
                                                                       9179, 20059, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 33109, 3, 6,
                                                                       31639, 19179, 31729, 9179,
                                                                       9239, 20159, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 33259, 3, 6,
                                                                       31729, 19239, 31819, 9239,
                                                                       9299, 20259, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 33409, 3, 6,
                                                                       31909, 19359, 32059, 9419,
                                                                       9509, 20359, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 33634, 3, 6,
                                                                       32059, 19459, 32209, 9509,
                                                                       9599, 20509, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 33859, 3, 6,
                                                                       32209, 19559, 32359, 9599,
                                                                       9689, 20659, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 34084, 3, 6,
                                                                       32359, 19659, 32509, 9689,
                                                                       9779, 20809, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 34309, 3, 6,
                                                                       32659, 19859, 32809, 9959,
                                                                       10049, 20959, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 34534, 3, 6,
                                                                       32809, 19959, 32959,
                                                                       10049, 10139, 21109,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 34759, 3, 6,
                                                                       32959, 20059, 33109,
                                                                       10139, 10229, 21259,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 34984, 3, 6,
                                                                       33109, 20159, 33259,
                                                                       10229, 10319, 21409,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35209, 0, 6,
                                                                       29959, 18059, 29974,
                                                                       10499, 10517, 21559,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35254, 0, 6,
                                                                       29974, 18069, 29989,
                                                                       10517, 10535, 21589,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35299, 0, 6,
                                                                       29989, 18079, 30004,
                                                                       10535, 10553, 21619,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35344, 0, 6,
                                                                       30004, 18089, 30019,
                                                                       10553, 10571, 21649,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35389, 0, 6,
                                                                       30019, 18099, 30034,
                                                                       10571, 10589, 21679,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35434, 0, 6,
                                                                       30034, 18109, 30049,
                                                                       10589, 10607, 21709,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35479, 0, 6,
                                                                       30049, 18119, 30064,
                                                                       10607, 10625, 21739,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35524, 0, 6,
                                                                       30079, 18139, 30094,
                                                                       10661, 10679, 21769,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35569, 0, 6,
                                                                       30094, 18149, 30109,
                                                                       10679, 10697, 21799,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35614, 0, 6,
                                                                       30109, 18159, 30124,
                                                                       10697, 10715, 21829,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35659, 0, 6,
                                                                       30124, 18169, 30139,
                                                                       10715, 10733, 21859,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35704, 0, 6,
                                                                       30139, 18179, 30154,
                                                                       10733, 10751, 21889,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35749, 0, 6,
                                                                       30154, 18189, 30169,
                                                                       10751, 10769, 21919,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 35794, 0, 6,
                                                                       30169, 18199, 30184,
                                                                       10769, 10787, 21949,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 35839, 0, 3, 6,
                                                                       30199, 18219, 30244,
                                                                       35209, 21559, 35254,
                                                                       10823, 10877, 21979,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 35974, 0, 3, 6,
                                                                       30244, 18249, 30289,
                                                                       35254, 21589, 35299,
                                                                       10877, 10931, 22069,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 36109, 0, 3, 6,
                                                                       30289, 18279, 30334,
                                                                       35299, 21619, 35344,
                                                                       10931, 10985, 22159,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 36244, 0, 3, 6,
                                                                       30334, 18309, 30379,
                                                                       35344, 21649, 35389,
                                                                       10985, 11039, 22249,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 36379, 0, 3, 6,
                                                                       30379, 18339, 30424,
                                                                       35389, 21679, 35434,
                                                                       11039, 11093, 22339,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 36514, 0, 3, 6,
                                                                       30424, 18369, 30469,
                                                                       35434, 21709, 35479,
                                                                       11093, 11147, 22429,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 36649, 0, 3, 6,
                                                                       30514, 18429, 30559,
                                                                       35524, 21769, 35569,
                                                                       11255, 11309, 22519,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 36784, 0, 3, 6,
                                                                       30559, 18459, 30604,
                                                                       35569, 21799, 35614,
                                                                       11309, 11363, 22609,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 36919, 0, 3, 6,
                                                                       30604, 18489, 30649,
                                                                       35614, 21829, 35659,
                                                                       11363, 11417, 22699,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 37054, 0, 3, 6,
                                                                       30649, 18519, 30694,
                                                                       35659, 21859, 35704,
                                                                       11417, 11471, 22789,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 37189, 0, 3, 6,
                                                                       30694, 18549, 30739,
                                                                       35704, 21889, 35749,
                                                                       11471, 11525, 22879,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 37324, 0, 3, 6,
                                                                       30739, 18579, 30784,
                                                                       35749, 21919, 35794,
                                                                       11525, 11579, 22969,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 37459, 0, 3, 6,
                                                                       30829, 18639, 30919,
                                                                       35839, 21979, 35974,
                                                                       11687, 11795, 23059,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 37729, 0, 3, 6,
                                                                       30919, 18699, 31009,
                                                                       35974, 22069, 36109,
                                                                       11795, 11903, 23239,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 37999, 0, 3, 6,
                                                                       31009, 18759, 31099,
                                                                       36109, 22159, 36244,
                                                                       11903, 12011, 23419,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 38269, 0, 3, 6,
                                                                       31099, 18819, 31189,
                                                                       36244, 22249, 36379,
                                                                       12011, 12119, 23599,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 38539, 0, 3, 6,
                                                                       31189, 18879, 31279,
                                                                       36379, 22339, 36514,
                                                                       12119, 12227, 23779,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 38809, 0, 3, 6,
                                                                       31369, 18999, 31459,
                                                                       36649, 22519, 36784,
                                                                       12443, 12551, 23959,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 39079, 0, 3, 6,
                                                                       31459, 19059, 31549,
                                                                       36784, 22609, 36919,
                                                                       12551, 12659, 24139,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 39349, 0, 3, 6,
                                                                       31549, 19119, 31639,
                                                                       36919, 22699, 37054,
                                                                       12659, 12767, 24319,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 39619, 0, 3, 6,
                                                                       31639, 19179, 31729,
                                                                       37054, 22789, 37189,
                                                                       12767, 12875, 24499,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 39889, 0, 3, 6,
                                                                       31729, 19239, 31819,
                                                                       37189, 22879, 37324,
                                                                       12875, 12983, 24679,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 40159, 0, 3, 6,
                                                                       31909, 19359, 32059,
                                                                       35839, 35974, 37459,
                                                                       23059, 37729, 13199,
                                                                       13379, 24859, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 40609, 0, 3, 6,
                                                                       32059, 19459, 32209,
                                                                       35974, 36109, 37729,
                                                                       23239, 37999, 13379,
                                                                       13559, 25159, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 41059, 0, 3, 6,
                                                                       32209, 19559, 32359,
                                                                       36109, 36244, 37999,
                                                                       23419, 38269, 13559,
                                                                       13739, 25459, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 41509, 0, 3, 6,
                                                                       32359, 19659, 32509,
                                                                       36244, 36379, 38269,
                                                                       23599, 38539, 13739,
                                                                       13919, 25759, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 41959, 0, 3, 6,
                                                                       32659, 19859, 32809,
                                                                       36649, 36784, 38809,
                                                                       23959, 39079, 14279,
                                                                       14459, 26059, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 42409, 0, 3, 6,
                                                                       32809, 19959, 32959,
                                                                       36784, 36919, 39079,
                                                                       24139, 39349, 14459,
                                                                       14639, 26359, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 42859, 0, 3, 6,
                                                                       32959, 20059, 33109,
                                                                       36919, 37054, 39349,
                                                                       24319, 39619, 14639,
                                                                       14819, 26659, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 43309, 0, 3, 6,
                                                                       33109, 20159, 33259,
                                                                       37054, 37189, 39619,
                                                                       24499, 39889, 14819,
                                                                       14999, 26959, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 43759, 0, 3, 6,
                                                                       33409, 20359, 33634,
                                                                       37459, 37729, 40159,
                                                                       24859, 40609, 15359,
                                                                       15629, 27259, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 44434, 0, 3, 6,
                                                                       33634, 20509, 33859,
                                                                       37729, 37999, 40609,
                                                                       25159, 41059, 15629,
                                                                       15899, 27709, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 45109, 0, 3, 6,
                                                                       33859, 20659, 34084,
                                                                       37999, 38269, 41059,
                                                                       25459, 41509, 15899,
                                                                       16169, 28159, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 45784, 0, 3, 6,
                                                                       34309, 20959, 34534,
                                                                       38809, 39079, 41959,
                                                                       26059, 42409, 16709,
                                                                       16979, 28609, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 46459, 0, 3, 6,
                                                                       34534, 21109, 34759,
                                                                       39079, 39349, 42409,
                                                                       26359, 42859, 16979,
                                                                       17249, 29059, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 47134, 0, 3, 6,
                                                                       34759, 21259, 34984,
                                                                       39349, 39619, 42859,
                                                                       26659, 43309, 17249,
                                                                       17519, 29509, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47809, 6, 18059,
                                                                       18069, 29989, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47830, 6, 18069,
                                                                       18079, 30004, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47851, 6, 18079,
                                                                       18089, 30019, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47872, 6, 18089,
                                                                       18099, 30034, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47893, 6, 18099,
                                                                       18109, 30049, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47914, 6, 18109,
                                                                       18119, 30064, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47935, 6, 18139,
                                                                       18149, 30109, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47956, 6, 18149,
                                                                       18159, 30124, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47977, 6, 18159,
                                                                       18169, 30139, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47998, 6, 18169,
                                                                       18179, 30154, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48019, 6, 18179,
                                                                       18189, 30169, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48040, 6, 18189,
                                                                       18199, 30184, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48061, 3, 6,
                                                                       47809, 29989, 47830,
                                                                       18219, 18249, 30289,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48124, 3, 6,
                                                                       47830, 30004, 47851,
                                                                       18249, 18279, 30334,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48187, 3, 6,
                                                                       47851, 30019, 47872,
                                                                       18279, 18309, 30379,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48250, 3, 6,
                                                                       47872, 30034, 47893,
                                                                       18309, 18339, 30424,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48313, 3, 6,
                                                                       47893, 30049, 47914,
                                                                       18339, 18369, 30469,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48376, 3, 6,
                                                                       47935, 30109, 47956,
                                                                       18429, 18459, 30604,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48439, 3, 6,
                                                                       47956, 30124, 47977,
                                                                       18459, 18489, 30649,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48502, 3, 6,
                                                                       47977, 30139, 47998,
                                                                       18489, 18519, 30694,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48565, 3, 6,
                                                                       47998, 30154, 48019,
                                                                       18519, 18549, 30739,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48628, 3, 6,
                                                                       48019, 30169, 48040,
                                                                       18549, 18579, 30784,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 48691, 3, 6,
                                                                       48061, 30289, 48124,
                                                                       18639, 18699, 31009,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 48817, 3, 6,
                                                                       48124, 30334, 48187,
                                                                       18699, 18759, 31099,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 48943, 3, 6,
                                                                       48187, 30379, 48250,
                                                                       18759, 18819, 31189,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 49069, 3, 6,
                                                                       48250, 30424, 48313,
                                                                       18819, 18879, 31279,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 49195, 3, 6,
                                                                       48376, 30604, 48439,
                                                                       18999, 19059, 31549,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 49321, 3, 6,
                                                                       48439, 30649, 48502,
                                                                       19059, 19119, 31639,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 49447, 3, 6,
                                                                       48502, 30694, 48565,
                                                                       19119, 19179, 31729,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 49573, 3, 6,
                                                                       48565, 30739, 48628,
                                                                       19179, 19239, 31819,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 49699, 3, 6,
                                                                       48691, 31009, 48817,
                                                                       19359, 19459, 32209,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 49909, 3, 6,
                                                                       48817, 31099, 48943,
                                                                       19459, 19559, 32359,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 50119, 3, 6,
                                                                       48943, 31189, 49069,
                                                                       19559, 19659, 32509,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 50329, 3, 6,
                                                                       49195, 31549, 49321,
                                                                       19859, 19959, 32959,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 50539, 3, 6,
                                                                       49321, 31639, 49447,
                                                                       19959, 20059, 33109,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 50749, 3, 6,
                                                                       49447, 31729, 49573,
                                                                       20059, 20159, 33259,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 50959, 3, 6,
                                                                       49699, 32209, 49909,
                                                                       20359, 20509, 33859,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 51274, 3, 6,
                                                                       49909, 32359, 50119,
                                                                       20509, 20659, 34084,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 51589, 3, 6,
                                                                       50329, 32959, 50539,
                                                                       20959, 21109, 34759,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 51904, 3, 6,
                                                                       50539, 33109, 50749,
                                                                       21109, 21259, 34984,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52219, 0, 6,
                                                                       47809, 29989, 47830,
                                                                       21559, 21589, 35299,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52282, 0, 6,
                                                                       47830, 30004, 47851,
                                                                       21589, 21619, 35344,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52345, 0, 6,
                                                                       47851, 30019, 47872,
                                                                       21619, 21649, 35389,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52408, 0, 6,
                                                                       47872, 30034, 47893,
                                                                       21649, 21679, 35434,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52471, 0, 6,
                                                                       47893, 30049, 47914,
                                                                       21679, 21709, 35479,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52534, 0, 6,
                                                                       47935, 30109, 47956,
                                                                       21769, 21799, 35614,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52597, 0, 6,
                                                                       47956, 30124, 47977,
                                                                       21799, 21829, 35659,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52660, 0, 6,
                                                                       47977, 30139, 47998,
                                                                       21829, 21859, 35704,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52723, 0, 6,
                                                                       47998, 30154, 48019,
                                                                       21859, 21889, 35749,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 52786, 0, 6,
                                                                       48019, 30169, 48040,
                                                                       21889, 21919, 35794,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 52849, 0, 3, 6,
                                                                       48061, 30289, 48124,
                                                                       52219, 35299, 52282,
                                                                       21979, 22069, 36109,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 53038, 0, 3, 6,
                                                                       48124, 30334, 48187,
                                                                       52282, 35344, 52345,
                                                                       22069, 22159, 36244,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 53227, 0, 3, 6,
                                                                       48187, 30379, 48250,
                                                                       52345, 35389, 52408,
                                                                       22159, 22249, 36379,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 53416, 0, 3, 6,
                                                                       48250, 30424, 48313,
                                                                       52408, 35434, 52471,
                                                                       22249, 22339, 36514,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 53605, 0, 3, 6,
                                                                       48376, 30604, 48439,
                                                                       52534, 35614, 52597,
                                                                       22519, 22609, 36919,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 53794, 0, 3, 6,
                                                                       48439, 30649, 48502,
                                                                       52597, 35659, 52660,
                                                                       22609, 22699, 37054,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 53983, 0, 3, 6,
                                                                       48502, 30694, 48565,
                                                                       52660, 35704, 52723,
                                                                       22699, 22789, 37189,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 54172, 0, 3, 6,
                                                                       48565, 30739, 48628,
                                                                       52723, 35749, 52786,
                                                                       22789, 22879, 37324,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 54361, 0, 3, 6,
                                                                       48691, 31009, 48817,
                                                                       52849, 36109, 53038,
                                                                       23059, 23239, 37999,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 54739, 0, 3, 6,
                                                                       48817, 31099, 48943,
                                                                       53038, 36244, 53227,
                                                                       23239, 23419, 38269,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 55117, 0, 3, 6,
                                                                       48943, 31189, 49069,
                                                                       53227, 36379, 53416,
                                                                       23419, 23599, 38539,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 55495, 0, 3, 6,
                                                                       49195, 31549, 49321,
                                                                       53605, 36919, 53794,
                                                                       23959, 24139, 39349,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 55873, 0, 3, 6,
                                                                       49321, 31639, 49447,
                                                                       53794, 37054, 53983,
                                                                       24139, 24319, 39619,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 56251, 0, 3, 6,
                                                                       49447, 31729, 49573,
                                                                       53983, 37189, 54172,
                                                                       24319, 24499, 39889,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 56629, 0, 3, 6,
                                                                       49699, 32209, 49909,
                                                                       52849, 53038, 54361,
                                                                       37999, 54739, 24859,
                                                                       25159, 41059, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 57259, 0, 3, 6,
                                                                       49909, 32359, 50119,
                                                                       53038, 53227, 54739,
                                                                       38269, 55117, 25159,
                                                                       25459, 41509, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 57889, 0, 3, 6,
                                                                       50329, 32959, 50539,
                                                                       53605, 53794, 55495,
                                                                       39349, 55873, 26059,
                                                                       26359, 42859, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 58519, 0, 3, 6,
                                                                       50539, 33109, 50749,
                                                                       53794, 53983, 55873,
                                                                       39619, 56251, 26359,
                                                                       26659, 43309, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgh_three_center_electron_repulsion_0(buffer, 59149, 0, 3, 6,
                                                                       50959, 33859, 51274,
                                                                       54361, 54739, 56629,
                                                                       41059, 57259, 27259,
                                                                       27709, 45109, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgh_three_center_electron_repulsion_0(buffer, 60094, 0, 3, 6,
                                                                       51589, 34759, 51904,
                                                                       55495, 55873, 57889,
                                                                       42859, 58519, 28609,
                                                                       29059, 47134, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61039, 6, 29959,
                                                                       29974, 47809, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61067, 6, 29974,
                                                                       29989, 47830, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61095, 6, 29989,
                                                                       30004, 47851, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61123, 6, 30004,
                                                                       30019, 47872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61151, 6, 30019,
                                                                       30034, 47893, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61179, 6, 30034,
                                                                       30049, 47914, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61207, 6, 30079,
                                                                       30094, 47935, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61235, 6, 30094,
                                                                       30109, 47956, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61263, 6, 30109,
                                                                       30124, 47977, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61291, 6, 30124,
                                                                       30139, 47998, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61319, 6, 30139,
                                                                       30154, 48019, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61347, 6, 30154,
                                                                       30169, 48040, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61375, 3, 6,
                                                                       61039, 47809, 61067,
                                                                       30199, 30244, 48061,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61459, 3, 6,
                                                                       61067, 47830, 61095,
                                                                       30244, 30289, 48124,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61543, 3, 6,
                                                                       61095, 47851, 61123,
                                                                       30289, 30334, 48187,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61627, 3, 6,
                                                                       61123, 47872, 61151,
                                                                       30334, 30379, 48250,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61711, 3, 6,
                                                                       61151, 47893, 61179,
                                                                       30379, 30424, 48313,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61795, 3, 6,
                                                                       61207, 47935, 61235,
                                                                       30514, 30559, 48376,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61879, 3, 6,
                                                                       61235, 47956, 61263,
                                                                       30559, 30604, 48439,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61963, 3, 6,
                                                                       61263, 47977, 61291,
                                                                       30604, 30649, 48502,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 62047, 3, 6,
                                                                       61291, 47998, 61319,
                                                                       30649, 30694, 48565,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 62131, 3, 6,
                                                                       61319, 48019, 61347,
                                                                       30694, 30739, 48628,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 62215, 3, 6,
                                                                       61375, 48061, 61459,
                                                                       30829, 30919, 48691,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 62383, 3, 6,
                                                                       61459, 48124, 61543,
                                                                       30919, 31009, 48817,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 62551, 3, 6,
                                                                       61543, 48187, 61627,
                                                                       31009, 31099, 48943,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 62719, 3, 6,
                                                                       61627, 48250, 61711,
                                                                       31099, 31189, 49069,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 62887, 3, 6,
                                                                       61795, 48376, 61879,
                                                                       31369, 31459, 49195,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 63055, 3, 6,
                                                                       61879, 48439, 61963,
                                                                       31459, 31549, 49321,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 63223, 3, 6,
                                                                       61963, 48502, 62047,
                                                                       31549, 31639, 49447,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 63391, 3, 6,
                                                                       62047, 48565, 62131,
                                                                       31639, 31729, 49573,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 63559, 3, 6,
                                                                       62215, 48691, 62383,
                                                                       31909, 32059, 49699,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 63839, 3, 6,
                                                                       62383, 48817, 62551,
                                                                       32059, 32209, 49909,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 64119, 3, 6,
                                                                       62551, 48943, 62719,
                                                                       32209, 32359, 50119,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 64399, 3, 6,
                                                                       62887, 49195, 63055,
                                                                       32659, 32809, 50329,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 64679, 3, 6,
                                                                       63055, 49321, 63223,
                                                                       32809, 32959, 50539,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 64959, 3, 6,
                                                                       63223, 49447, 63391,
                                                                       32959, 33109, 50749,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 65239, 3, 6,
                                                                       63559, 49699, 63839,
                                                                       33409, 33634, 50959,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 65659, 3, 6,
                                                                       63839, 49909, 64119,
                                                                       33634, 33859, 51274,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 66079, 3, 6,
                                                                       64399, 50329, 64679,
                                                                       34309, 34534, 51589,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 66499, 3, 6,
                                                                       64679, 50539, 64959,
                                                                       34534, 34759, 51904,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 66919, 0, 6,
                                                                       61039, 47809, 61067,
                                                                       35209, 35254, 52219,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 67003, 0, 6,
                                                                       61067, 47830, 61095,
                                                                       35254, 35299, 52282,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 67087, 0, 6,
                                                                       61095, 47851, 61123,
                                                                       35299, 35344, 52345,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 67171, 0, 6,
                                                                       61123, 47872, 61151,
                                                                       35344, 35389, 52408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 67255, 0, 6,
                                                                       61151, 47893, 61179,
                                                                       35389, 35434, 52471,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 67339, 0, 6,
                                                                       61207, 47935, 61235,
                                                                       35524, 35569, 52534,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 67423, 0, 6,
                                                                       61235, 47956, 61263,
                                                                       35569, 35614, 52597,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 67507, 0, 6,
                                                                       61263, 47977, 61291,
                                                                       35614, 35659, 52660,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 67591, 0, 6,
                                                                       61291, 47998, 61319,
                                                                       35659, 35704, 52723,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 67675, 0, 6,
                                                                       61319, 48019, 61347,
                                                                       35704, 35749, 52786,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 67759, 0, 3, 6,
                                                                       61375, 48061, 61459,
                                                                       66919, 52219, 67003,
                                                                       35839, 35974, 52849,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 68011, 0, 3, 6,
                                                                       61459, 48124, 61543,
                                                                       67003, 52282, 67087,
                                                                       35974, 36109, 53038,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 68263, 0, 3, 6,
                                                                       61543, 48187, 61627,
                                                                       67087, 52345, 67171,
                                                                       36109, 36244, 53227,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 68515, 0, 3, 6,
                                                                       61627, 48250, 61711,
                                                                       67171, 52408, 67255,
                                                                       36244, 36379, 53416,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 68767, 0, 3, 6,
                                                                       61795, 48376, 61879,
                                                                       67339, 52534, 67423,
                                                                       36649, 36784, 53605,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 69019, 0, 3, 6,
                                                                       61879, 48439, 61963,
                                                                       67423, 52597, 67507,
                                                                       36784, 36919, 53794,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 69271, 0, 3, 6,
                                                                       61963, 48502, 62047,
                                                                       67507, 52660, 67591,
                                                                       36919, 37054, 53983,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 69523, 0, 3, 6,
                                                                       62047, 48565, 62131,
                                                                       67591, 52723, 67675,
                                                                       37054, 37189, 54172,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 69775, 0, 3, 6,
                                                                       62215, 48691, 62383,
                                                                       67759, 52849, 68011,
                                                                       37459, 37729, 54361,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 70279, 0, 3, 6,
                                                                       62383, 48817, 62551,
                                                                       68011, 53038, 68263,
                                                                       37729, 37999, 54739,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 70783, 0, 3, 6,
                                                                       62551, 48943, 62719,
                                                                       68263, 53227, 68515,
                                                                       37999, 38269, 55117,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 71287, 0, 3, 6,
                                                                       62887, 49195, 63055,
                                                                       68767, 53605, 69019,
                                                                       38809, 39079, 55495,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 71791, 0, 3, 6,
                                                                       63055, 49321, 63223,
                                                                       69019, 53794, 69271,
                                                                       39079, 39349, 55873,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 72295, 0, 3, 6,
                                                                       63223, 49447, 63391,
                                                                       69271, 53983, 69523,
                                                                       39349, 39619, 56251,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfi_three_center_electron_repulsion_0(buffer, 72799, 0, 3, 6,
                                                                       63559, 49699, 63839,
                                                                       67759, 68011, 69775,
                                                                       54361, 70279, 40159,
                                                                       40609, 56629, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfi_three_center_electron_repulsion_0(buffer, 73639, 0, 3, 6,
                                                                       63839, 49909, 64119,
                                                                       68011, 68263, 70279,
                                                                       54739, 70783, 40609,
                                                                       41059, 57259, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfi_three_center_electron_repulsion_0(buffer, 74479, 0, 3, 6,
                                                                       64399, 50329, 64679,
                                                                       68767, 69019, 71287,
                                                                       55495, 71791, 41959,
                                                                       42409, 57889, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfi_three_center_electron_repulsion_0(buffer, 75319, 0, 3, 6,
                                                                       64679, 50539, 64959,
                                                                       69019, 69271, 71791,
                                                                       55873, 72295, 42409,
                                                                       42859, 58519, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgi_three_center_electron_repulsion_0(buffer, 76159, 0, 3, 6,
                                                                       65239, 50959, 65659,
                                                                       69775, 70279, 72799,
                                                                       56629, 73639, 43759,
                                                                       44434, 59149, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgi_three_center_electron_repulsion_0(buffer, 77419, 0, 3, 6,
                                                                       66079, 51589, 66499,
                                                                       71287, 71791, 74479,
                                                                       57889, 75319, 45784,
                                                                       46459, 60094, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 78679, 77419, 1, 420, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 79099, 77419, 1, 420, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 79519, 77419, 1, 420, ncols, alpha);

                    simdgeo::geom_s_x(buffer, 79939, 76159, 1, 420, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 80359, 76159, 1, 420, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 80779, 76159, 1, 420, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 81199, 78679, 2520, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 83719, 81199, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 83719, 13, nmax);

        simdtrf::transform_i_inner(buffer, 83719, 81619, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 117 * nvalues + n * npairs, nvalues, buffer, 83719,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 83719, 82039, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 234 * nvalues + n * npairs, nvalues, buffer, 83719,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 83719, 82459, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 351 * nvalues + n * npairs, nvalues, buffer, 83719,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 83719, 82879, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 468 * nvalues + n * npairs, nvalues, buffer, 83719,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 83719, 83299, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 585 * nvalues + n * npairs, nvalues, buffer, 83719,
                                   13, nmax);
    }

    for (size_t m = 0; m < 702; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
