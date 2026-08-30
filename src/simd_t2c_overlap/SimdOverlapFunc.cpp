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


#include "SimdOverlapFunc.hpp"

#include <string>

#include "ErrorHandler.hpp"
#include "SimdOverlapRecDD.hpp"
#include "SimdOverlapRecDF.hpp"
#include "SimdOverlapRecDG.hpp"
#include "SimdOverlapRecDH.hpp"
#include "SimdOverlapRecDI.hpp"
#include "SimdOverlapRecDP.hpp"
#include "SimdOverlapRecDS.hpp"
#include "SimdOverlapRecFD.hpp"
#include "SimdOverlapRecFF.hpp"
#include "SimdOverlapRecFG.hpp"
#include "SimdOverlapRecFH.hpp"
#include "SimdOverlapRecFI.hpp"
#include "SimdOverlapRecFP.hpp"
#include "SimdOverlapRecFS.hpp"
#include "SimdOverlapRecGD.hpp"
#include "SimdOverlapRecGF.hpp"
#include "SimdOverlapRecGG.hpp"
#include "SimdOverlapRecGH.hpp"
#include "SimdOverlapRecGI.hpp"
#include "SimdOverlapRecGP.hpp"
#include "SimdOverlapRecGS.hpp"
#include "SimdOverlapRecHD.hpp"
#include "SimdOverlapRecHF.hpp"
#include "SimdOverlapRecHG.hpp"
#include "SimdOverlapRecHH.hpp"
#include "SimdOverlapRecHI.hpp"
#include "SimdOverlapRecHP.hpp"
#include "SimdOverlapRecHS.hpp"
#include "SimdOverlapRecID.hpp"
#include "SimdOverlapRecIF.hpp"
#include "SimdOverlapRecIG.hpp"
#include "SimdOverlapRecIH.hpp"
#include "SimdOverlapRecII.hpp"
#include "SimdOverlapRecIP.hpp"
#include "SimdOverlapRecIS.hpp"
#include "SimdOverlapRecPD.hpp"
#include "SimdOverlapRecPF.hpp"
#include "SimdOverlapRecPG.hpp"
#include "SimdOverlapRecPH.hpp"
#include "SimdOverlapRecPI.hpp"
#include "SimdOverlapRecPP.hpp"
#include "SimdOverlapRecPS.hpp"
#include "SimdOverlapRecSD.hpp"
#include "SimdOverlapRecSF.hpp"
#include "SimdOverlapRecSG.hpp"
#include "SimdOverlapRecSH.hpp"
#include "SimdOverlapRecSI.hpp"
#include "SimdOverlapRecSP.hpp"
#include "SimdOverlapRecSS.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_overlap(double                         *values,
                const size_t                    nvalues,
                const CBasisFunction           &bra,
                const CBasisFunction           &ket,
                const std::vector<CSimdMatrix> &harmonics,
                const CSimdMatrix              &coordinates,
                const double                    threshold) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    // NOTE: the kernel of the overlap of two S type functions needs no solid
    // harmonics, as the harmonics of angular momentum zero are one for every
    // atom pair.

    if ((lbra == 0) && (lket == 0))
    {
        compute_ss_overlap(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    // NOTE: the kernels of the combinations with one S type function take the
    // solid harmonics of the angular momentum of the other side alone, as their
    // integrals are a single term.

    if ((lbra == 0) && (lket == 1))
    {
        compute_sp_overlap(values, nvalues, bra, ket, harmonics[0], coordinates, threshold);

        return;
    }

    if ((lbra == 0) && (lket == 2))
    {
        compute_sd_overlap(values, nvalues, bra, ket, harmonics[1], coordinates, threshold);

        return;
    }

    if ((lbra == 0) && (lket == 3))
    {
        compute_sf_overlap(values, nvalues, bra, ket, harmonics[2], coordinates, threshold);

        return;
    }

    if ((lbra == 0) && (lket == 4))
    {
        compute_sg_overlap(values, nvalues, bra, ket, harmonics[3], coordinates, threshold);

        return;
    }

    if ((lbra == 0) && (lket == 5))
    {
        compute_sh_overlap(values, nvalues, bra, ket, harmonics[4], coordinates, threshold);

        return;
    }

    if ((lbra == 0) && (lket == 6))
    {
        compute_si_overlap(values, nvalues, bra, ket, harmonics[5], coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 0))
    {
        compute_ps_overlap(values, nvalues, bra, ket, harmonics[0], coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 0))
    {
        compute_ds_overlap(values, nvalues, bra, ket, harmonics[1], coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 0))
    {
        compute_fs_overlap(values, nvalues, bra, ket, harmonics[2], coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 0))
    {
        compute_gs_overlap(values, nvalues, bra, ket, harmonics[3], coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 0))
    {
        compute_hs_overlap(values, nvalues, bra, ket, harmonics[4], coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 0))
    {
        compute_is_overlap(values, nvalues, bra, ket, harmonics[5], coordinates, threshold);

        return;
    }

    // NOTE: the remaining kernels take all the solid harmonics, as their terms
    // reach every angular momentum up to the sum of the two.

    if ((lbra == 1) && (lket == 1))
    {
        compute_pp_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 2))
    {
        compute_pd_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 3))
    {
        compute_pf_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 4))
    {
        compute_pg_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 5))
    {
        compute_ph_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 1) && (lket == 6))
    {
        compute_pi_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 1))
    {
        compute_dp_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 2))
    {
        compute_dd_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 3))
    {
        compute_df_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 4))
    {
        compute_dg_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 5))
    {
        compute_dh_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 2) && (lket == 6))
    {
        compute_di_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 1))
    {
        compute_fp_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 2))
    {
        compute_fd_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 3))
    {
        compute_ff_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 4))
    {
        compute_fg_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 5))
    {
        compute_fh_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 3) && (lket == 6))
    {
        compute_fi_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 1))
    {
        compute_gp_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 2))
    {
        compute_gd_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 3))
    {
        compute_gf_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 4))
    {
        compute_gg_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 5))
    {
        compute_gh_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 4) && (lket == 6))
    {
        compute_gi_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 1))
    {
        compute_hp_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 2))
    {
        compute_hd_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 3))
    {
        compute_hf_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 4))
    {
        compute_hg_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 5))
    {
        compute_hh_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 5) && (lket == 6))
    {
        compute_hi_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 1))
    {
        compute_ip_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 2))
    {
        compute_id_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 3))
    {
        compute_if_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 4))
    {
        compute_ig_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 5))
    {
        compute_ih_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    if ((lbra == 6) && (lket == 6))
    {
        compute_ii_overlap(values, nvalues, bra, ket, harmonics, coordinates, threshold);

        return;
    }

    errors::assertMsgCritical(
        false, std::string("SimdOverlapFunc.compute_overlap: Overlap integrals of the requested angular momenta are not implemented"));
}

}  // namespace simdovl
