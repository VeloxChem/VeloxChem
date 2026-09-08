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
compute_overlap(double               *values,
                const size_t          nvalues,
                const CBasisFunction &bra,
                const CBasisFunction &ket,
                const CSimdMatrix    &coordinates,
                const double          threshold) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    // NOTE: the two orders of a combination are separate kernels, as the recurrence
    // builds the angular momentum on ket side and transfers it to bra side, so the
    // work differs with the order even where the integrals do not.

    // NOTE: the combination is dispatched on the pair of angular momenta taken as a
    // single index, so the compiler forms one jump table rather than a chain of
    // comparisons which the combinations of high angular momentum reach last.

    switch (lbra * 7 + lket)
    {
        case  0:
            compute_ss_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  1:
            compute_sp_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  2:
            compute_sd_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  3:
            compute_sf_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  4:
            compute_sg_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  5:
            compute_sh_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  6:
            compute_si_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  7:
            compute_ps_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  8:
            compute_pp_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case  9:
            compute_pd_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 10:
            compute_pf_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 11:
            compute_pg_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 12:
            compute_ph_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 13:
            compute_pi_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 14:
            compute_ds_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 15:
            compute_dp_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 16:
            compute_dd_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 17:
            compute_df_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 18:
            compute_dg_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 19:
            compute_dh_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 20:
            compute_di_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 21:
            compute_fs_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 22:
            compute_fp_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 23:
            compute_fd_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 24:
            compute_ff_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 25:
            compute_fg_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 26:
            compute_fh_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 27:
            compute_fi_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 28:
            compute_gs_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 29:
            compute_gp_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 30:
            compute_gd_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 31:
            compute_gf_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 32:
            compute_gg_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 33:
            compute_gh_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 34:
            compute_gi_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 35:
            compute_hs_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 36:
            compute_hp_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 37:
            compute_hd_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 38:
            compute_hf_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 39:
            compute_hg_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 40:
            compute_hh_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 41:
            compute_hi_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 42:
            compute_is_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 43:
            compute_ip_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 44:
            compute_id_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 45:
            compute_if_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 46:
            compute_ig_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 47:
            compute_ih_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        case 48:
            compute_ii_overlap(values, nvalues, bra, ket, coordinates, threshold);
            return;

        default:
            break;
    }

    // NOTE: the kernels are generated up to angular momentum six, so a combination
    // above it stops rather than leaving the values of the sparsity pattern
    // unwritten, which is what a caller would otherwise read as integrals.

    errors::assertMsgCritical(false, std::string("SimdOverlapFunc.compute_overlap: Overlap integrals are not implemented"));
}

}  // namespace simdovl
