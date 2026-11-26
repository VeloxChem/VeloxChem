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


#include "SimdThreeCenterElectronRepulsionGeom010RecGDD.hpp"

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
#include "SimdGeometryI1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDD.hpp"
#include "SimdTransferDF.hpp"
#include "SimdTransferFD.hpp"
#include "SimdTransferGeom010XDD.hpp"
#include "SimdTransferGeom010XDF.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XFD.hpp"
#include "SimdTransferGeom010XFF.hpp"
#include "SimdTransferGeom010XGD.hpp"
#include "SimdTransferGeom010XPD.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010YDD.hpp"
#include "SimdTransferGeom010YDF.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YFD.hpp"
#include "SimdTransferGeom010YFF.hpp"
#include "SimdTransferGeom010YGD.hpp"
#include "SimdTransferGeom010YPD.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010ZDD.hpp"
#include "SimdTransferGeom010ZDF.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZFD.hpp"
#include "SimdTransferGeom010ZFF.hpp"
#include "SimdTransferGeom010ZGD.hpp"
#include "SimdTransferGeom010ZPD.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_gdd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_gdd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 18699, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 675 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 18699, 4987, 3072, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto beta = b_exps[j];

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pb(buffer, coordinates, 0, nmax, fb);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 9, ncols,
                                                             fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 17, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 20, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 44, 0, 3, 7, 8,
                                                                       17, 20, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 50, 0, 3, 8, 9,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 56, 0, 3, 9, 10,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 62, 0, 3, 10, 11,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 68, 0, 3, 11, 12,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 74, 0, 3, 12, 13,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 80, 0, 3, 13, 14,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 86, 0, 3, 14, 15,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 92, 0, 3, 17, 20,
                                                                       44, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 102, 0, 3, 20, 23,
                                                                       50, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 112, 0, 3, 23, 26,
                                                                       56, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 122, 0, 3, 26, 29,
                                                                       62, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 132, 0, 3, 29, 32,
                                                                       68, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 142, 0, 3, 32, 35,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 152, 0, 3, 35, 38,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 162, 0, 3, 44, 50,
                                                                       92, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 177, 0, 3, 50, 56,
                                                                       102, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 192, 0, 3, 56, 62,
                                                                       112, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 207, 0, 3, 62, 68,
                                                                       122, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 222, 0, 3, 68, 74,
                                                                       132, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 237, 0, 3, 74, 80,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 252, 0, 3, 92,
                                                                       102, 162, 177, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 102,
                                                                       112, 177, 192, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 294, 0, 3, 112,
                                                                       122, 192, 207, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 315, 0, 3, 122,
                                                                       132, 207, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 336, 0, 3, 132,
                                                                       142, 222, 237, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 357, 0, 3, 162,
                                                                       177, 252, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 385, 0, 3, 177,
                                                                       192, 273, 294, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 413, 0, 3, 192,
                                                                       207, 294, 315, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 441, 0, 3, 207,
                                                                       222, 315, 336, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 469, 0, 3, 252,
                                                                       273, 357, 385, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 505, 0, 3, 273,
                                                                       294, 385, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 541, 0, 3, 294,
                                                                       315, 413, 441, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 577, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 580, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 583, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 586, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 589, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 592, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 595, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 598, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 601, 3, 9, 23,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 610, 3, 10, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 619, 3, 11, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 628, 3, 12, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 637, 3, 13, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 646, 3, 14, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 655, 3, 15, 41,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 664, 3, 23, 56,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 682, 3, 26, 62,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 700, 3, 29, 68,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 718, 3, 32, 74,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 736, 3, 35, 80,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 754, 3, 38, 86,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 772, 3, 56, 112,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 802, 3, 62, 122,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 832, 3, 68, 132,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 862, 3, 74, 142,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 892, 3, 80, 152,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 922, 3, 112, 192,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 967, 3, 122, 207,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1012, 3, 132, 222,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1057, 3, 142, 237,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1102, 3, 192, 294,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1165, 3, 207, 315,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 1228, 3, 222, 336,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1291, 3, 294, 413,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 1375, 3, 315, 441,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 1459, 3, 413, 541,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1567, 3, 7, 8,
                                                                       577, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1573, 3, 8, 9,
                                                                       580, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1579, 3, 9, 10,
                                                                       583, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1585, 3, 10, 11,
                                                                       586, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1591, 3, 11, 12,
                                                                       589, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1597, 3, 12, 13,
                                                                       592, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1603, 3, 13, 14,
                                                                       595, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1609, 3, 14, 15,
                                                                       598, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1615, 0, 3, 1567,
                                                                       577, 1573, 601, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1633, 0, 3, 1573,
                                                                       580, 1579, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1651, 0, 3, 1579,
                                                                       583, 1585, 619, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1669, 0, 3, 1585,
                                                                       586, 1591, 628, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1687, 0, 3, 1591,
                                                                       589, 1597, 637, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1705, 0, 3, 1597,
                                                                       592, 1603, 646, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1723, 0, 3, 1603,
                                                                       595, 1609, 655, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1741, 0, 3, 1615,
                                                                       601, 1633, 44, 50, 664,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1777, 0, 3, 1633,
                                                                       610, 1651, 50, 56, 682,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1813, 0, 3, 1651,
                                                                       619, 1669, 56, 62, 700,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1849, 0, 3, 1669,
                                                                       628, 1687, 62, 68, 718,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1885, 0, 3, 1687,
                                                                       637, 1705, 68, 74, 736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1921, 0, 3, 1705,
                                                                       646, 1723, 74, 80, 754,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 1957, 0, 3, 1741,
                                                                       664, 1777, 92, 102, 772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2017, 0, 3, 1777,
                                                                       682, 1813, 102, 112, 802,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2077, 0, 3, 1813,
                                                                       700, 1849, 112, 122, 832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2137, 0, 3, 1849,
                                                                       718, 1885, 122, 132, 862,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2197, 0, 3, 1885,
                                                                       736, 1921, 132, 142, 892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2257, 0, 3, 1957,
                                                                       772, 2017, 162, 177, 922,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2347, 0, 3, 2017,
                                                                       802, 2077, 177, 192, 967,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2437, 0, 3, 2077,
                                                                       832, 2137, 192, 207, 1012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 2527, 0, 3, 2137,
                                                                       862, 2197, 207, 222, 1057,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2617, 0, 3, 2257,
                                                                       922, 2347, 252, 273, 1102,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2743, 0, 3, 2347,
                                                                       967, 2437, 273, 294, 1165,
                                                                       ncols, gamma, p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 2869, 0, 3, 2437,
                                                                       1012, 2527, 294, 315,
                                                                       1228, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 2995, 0, 3, 2617,
                                                                       1102, 2743, 357, 385,
                                                                       1291, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 3163, 0, 3, 2743,
                                                                       1165, 2869, 385, 413,
                                                                       1375, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 3331, 0, 3, 2995,
                                                                       1291, 3163, 469, 505,
                                                                       1459, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_d_x(buffer, 3547, 1615, 1957, 1, 6, ncols, beta);

                    simdgeo::geom_d_y(buffer, 3583, 1615, 1957, 1, 6, ncols, beta);

                    simdgeo::geom_d_z(buffer, 3619, 1615, 1957, 1, 6, ncols, beta);

                    simdgeo::geom_f_x(buffer, 3655, 1741, 2257, 1, 6, ncols, beta);

                    simdgeo::geom_f_y(buffer, 3715, 1741, 2257, 1, 6, ncols, beta);

                    simdgeo::geom_f_z(buffer, 3775, 1741, 2257, 1, 6, ncols, beta);

                    simdgeo::geom_g_x(buffer, 3835, 1957, 2617, 1, 6, ncols, beta);

                    simdgeo::geom_g_y(buffer, 3925, 1957, 2617, 1, 6, ncols, beta);

                    simdgeo::geom_g_z(buffer, 4015, 1957, 2617, 1, 6, ncols, beta);

                    simdgeo::geom_h_x(buffer, 4105, 2257, 2995, 1, 6, ncols, beta);

                    simdgeo::geom_h_y(buffer, 4231, 2257, 2995, 1, 6, ncols, beta);

                    simdgeo::geom_h_z(buffer, 4357, 2257, 2995, 1, 6, ncols, beta);

                    simdgeo::geom_i_x(buffer, 4483, 2617, 3331, 1, 6, ncols, beta);

                    simdgeo::geom_i_y(buffer, 4651, 2617, 3331, 1, 6, ncols, beta);

                    simdgeo::geom_i_z(buffer, 4819, 2617, 3331, 1, 6, ncols, beta);

                    simdfunc::contract_primitives(buffer, 4987, 3547, 36, ncols);

                    simdfunc::contract_primitives(buffer, 5053, 3583, 36, ncols);

                    simdfunc::contract_primitives(buffer, 5119, 3619, 36, ncols);

                    simdfunc::contract_primitives(buffer, 5185, 1741, 36, ncols);

                    simdfunc::contract_primitives(buffer, 5251, 3655, 60, ncols);

                    simdfunc::contract_primitives(buffer, 5361, 3715, 60, ncols);

                    simdfunc::contract_primitives(buffer, 5471, 3775, 60, ncols);

                    simdfunc::contract_primitives(buffer, 5581, 1957, 60, ncols);

                    simdfunc::contract_primitives(buffer, 5691, 3835, 90, ncols);

                    simdfunc::contract_primitives(buffer, 5856, 3925, 90, ncols);

                    simdfunc::contract_primitives(buffer, 6021, 4015, 90, ncols);

                    simdfunc::contract_primitives(buffer, 6186, 2257, 90, ncols);

                    simdfunc::contract_primitives(buffer, 6351, 4105, 126, ncols);

                    simdfunc::contract_primitives(buffer, 6582, 4231, 126, ncols);

                    simdfunc::contract_primitives(buffer, 6813, 4357, 126, ncols);

                    simdfunc::contract_primitives(buffer, 7044, 2617, 126, ncols);

                    simdfunc::contract_primitives(buffer, 7275, 4483, 168, ncols);

                    simdfunc::contract_primitives(buffer, 7583, 4651, 168, ncols);

                    simdfunc::contract_primitives(buffer, 7891, 4819, 168, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 5023, 4987, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5089, 5053, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5155, 5119, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5221, 5185, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5311, 5251, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5421, 5361, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5531, 5471, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5641, 5581, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5781, 5691, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5946, 5856, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 6111, 6021, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 6276, 6186, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 6477, 6351, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 6708, 6582, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 6939, 6813, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 7170, 7044, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 7443, 7275, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 7751, 7583, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 8059, 7891, 28, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 8199, 5023, 5221, 5311, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 8289, 5089, 5221, 5421, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 8379, 5155, 5221, 5531, 5, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 8469, 5221, 5641, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 8559, 5311, 5641, 5781, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 8709, 5421, 5641, 5946, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 8859, 5531, 5641, 6111, 5, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 9009, 5641, 6276, 5, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 9159, 5781, 6276, 6477, 5, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 9384, 5946, 6276, 6708, 5, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 9609, 6111, 6276, 6939, 5, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 9834, 6276, 7170, 5, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 10059, 6477, 7170, 7443, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 10374, 6708, 7170, 7751, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 10689, 6939, 7170, 8059, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 11004, 8199, 8469, 8559, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 11184, 8289, 8469, 8709, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 11364, 8379, 8469, 8859, 5,
                                          nmax);

        simdtrf::compute_hrr_dd(buffer, coordinates, 11544, 8469, 9009, 5, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 11724, 8559, 9009, 9159, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 12024, 8709, 9009, 9384, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 12324, 8859, 9009, 9609, 5,
                                          nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 12624, 9009, 9834, 5, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 12924, 9159, 9834, 10059, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 13374, 9384, 9834, 10374, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 13824, 9609, 9834, 10689, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fd_out_of_second(buffer, coordinates, 14274, 11004, 11544,
                                                        11724, 5, nmax);

        simdtrf::compute_hrr_geom_010y_fd_out_of_second(buffer, coordinates, 14574, 11184, 11544,
                                                        12024, 5, nmax);

        simdtrf::compute_hrr_geom_010z_fd_out_of_second(buffer, coordinates, 14874, 11364, 11544,
                                                        12324, 5, nmax);

        simdtrf::compute_hrr_fd_out_of_second(buffer, coordinates, 15174, 11544, 12624, 5,
                                              nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 15474, 11724, 12624, 12924, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 15974, 12024, 12624, 13374, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 16474, 12324, 12624, 13824, 5,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gd_out_of_second(buffer, coordinates, 16974, 14274, 15174,
                                                        15474, 5, nmax);

        simdtrf::compute_hrr_geom_010y_gd_out_of_second(buffer, coordinates, 17424, 14574, 15174,
                                                        15974, 5, nmax);

        simdtrf::compute_hrr_geom_010z_gd_out_of_second(buffer, coordinates, 17874, 14874, 15174,
                                                        16474, 5, nmax);

        simdtrf::transform_d_inner(buffer, 18324, 16974, 15, 5, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 18324, 25, nmax);

        simdtrf::transform_d_inner(buffer, 18324, 17424, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 225 * nvalues + n * npairs, nvalues, buffer, 18324,
                                   25, nmax);

        simdtrf::transform_d_inner(buffer, 18324, 17874, 15, 5, nmax);

        simdtrf::transform_g_outer(values + 450 * nvalues + n * npairs, nvalues, buffer, 18324,
                                   25, nmax);
    }

    for (size_t m = 0; m < 675; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
