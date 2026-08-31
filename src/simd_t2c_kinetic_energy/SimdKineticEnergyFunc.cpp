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


#include "SimdKineticEnergyFunc.hpp"

#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"

namespace simdkin {  // simdkin namespace

auto
one_center_kinetic_energy(const CBasisFunction &bra, const CBasisFunction &ket) -> double
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    if (lbra != lket) return 0.0;

    // NOTE: the kinetic energy of two primitives on the same center is their
    // overlap times two l plus three times their reduced exponent, which is the
    // bracket of the integrals of separated atoms at zero separation.

    const auto ffact = static_cast<double>(2 * lbra + 3);

    const auto fdfact = mathfunc::double_factorial(2 * lbra - 1);

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    auto fsum = 0.0;

    for (size_t i = 0; i < a_exps.size(); i++)
    {
        for (size_t j = 0; j < b_exps.size(); j++)
        {
            const auto fab = 1.0 / (a_exps[i] + b_exps[j]);

            const auto fmu = a_exps[i] * b_exps[j] * fab;

            auto fovl = mathconst::pi_value() * fab;

            fovl = a_norms[i] * b_norms[j] * fovl * std::sqrt(fovl);

            // NOTE: repeated multiplication is used, as the angular momentum is
            // small.

            const auto fhalf = 0.5 * fab;

            for (int k = 0; k < lbra; k++)
            {
                fovl *= fhalf;
            }

            fsum += ffact * fmu * fdfact * fovl;
        }
    }

    return fsum;
}

auto
compute_kinetic_energy(double                         *values,
                       const size_t                    nvalues,
                       const CBasisFunction           &bra,
                       const CBasisFunction           &ket,
                       const std::vector<CSimdMatrix> &harmonics,
                       const CSimdMatrix              &coordinates,
                       const double                    threshold) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    // NOTE: the kernel of two S type functions needs no solid harmonics, as the
    // harmonics of angular momentum zero are one for every atom pair.

    if ((lbra == 0) && (lket == 0))
    {
        compute_ss_kinetic_energy(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    // NOTE: the kernels of the combinations with one S type function take the
    // solid harmonics of the angular momentum of the other side alone, as their
    // integrals are a single term.

    if ((lbra == 0) && (lket == 1))
    {
        compute_sp_kinetic_energy(values, nvalues, bra, ket, harmonics[0], coordinates, threshold);

        return;
    }

    if ((lbra == 0) && (lket == 2))
    {
        compute_sd_kinetic_energy(values, nvalues, bra, ket, harmonics[1], coordinates, threshold);

        return;
    }

    if ((lbra == 0) && (lket == 3))
    {
        compute_sf_kinetic_energy(values, nvalues, bra, ket, harmonics[2], coordinates, threshold);

        return;
    }

    if ((lbra == 0) && (lket == 4))
    {
        compute_sg_kinetic_energy(values, nvalues, bra, ket, harmonics[3], coordinates, threshold);

        return;
    }

    if ((lbra == 0) && (lket == 5))
    {
        compute_sh_kinetic_energy(values, nvalues, bra, ket, harmonics[4], coordinates, threshold);

        return;
    }

    if ((lbra == 0) && (lket == 6))
    {
        compute_si_kinetic_energy(values, nvalues, bra, ket, harmonics[5], coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 0))
    {
        compute_ps_kinetic_energy(values, nvalues, bra, ket, harmonics[0], coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 0))
    {
        compute_ds_kinetic_energy(values, nvalues, bra, ket, harmonics[1], coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 0))
    {
        compute_fs_kinetic_energy(values, nvalues, bra, ket, harmonics[2], coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 0))
    {
        compute_gs_kinetic_energy(values, nvalues, bra, ket, harmonics[3], coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 0))
    {
        compute_hs_kinetic_energy(values, nvalues, bra, ket, harmonics[4], coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 0))
    {
        compute_is_kinetic_energy(values, nvalues, bra, ket, harmonics[5], coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 1))
    {
        compute_pp_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 2))
    {
        compute_pd_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 3))
    {
        compute_pf_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 4))
    {
        compute_pg_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 5))
    {
        compute_ph_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 6))
    {
        compute_pi_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 1))
    {
        compute_dp_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 2))
    {
        compute_dd_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 3))
    {
        compute_df_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 4))
    {
        compute_dg_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 5))
    {
        compute_dh_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 6))
    {
        compute_di_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 1))
    {
        compute_fp_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 2))
    {
        compute_fd_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 3))
    {
        compute_ff_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 4))
    {
        compute_fg_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 5))
    {
        compute_fh_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 6))
    {
        compute_fi_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 1))
    {
        compute_gp_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 2))
    {
        compute_gd_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 3))
    {
        compute_gf_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 4))
    {
        compute_gg_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 5))
    {
        compute_gh_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 6))
    {
        compute_gi_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 1))
    {
        compute_hp_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 2))
    {
        compute_hd_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 3))
    {
        compute_hf_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 4))
    {
        compute_hg_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 5))
    {
        compute_hh_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 6))
    {
        compute_hi_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 1))
    {
        compute_ip_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 2))
    {
        compute_id_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 3))
    {
        compute_if_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 4))
    {
        compute_ig_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 5))
    {
        compute_ih_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 6))
    {
        compute_ii_kinetic_energy(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    errors::assertMsgCritical(false,
                              std::string("SimdKineticEnergyFunc.compute_kinetic_energy: Combination of angular momenta is not supported"));
}

}  // namespace simdkin
