#ifndef ThreeCenterOverlapFunc_hpp
#define ThreeCenterOverlapFunc_hpp

#include <array>

#include "GtoBlock.hpp"
#include "ThreeCenterOverlapSumRecGG.hpp"
#include "ThreeCenterOverlapSumRecGF.hpp"
#include "ThreeCenterOverlapSumRecGD.hpp"
#include "ThreeCenterOverlapSumRecGP.hpp"
#include "ThreeCenterOverlapSumRecGS.hpp"
#include "ThreeCenterOverlapSumRecFG.hpp"
#include "ThreeCenterOverlapSumRecFF.hpp"
#include "ThreeCenterOverlapSumRecFD.hpp"
#include "ThreeCenterOverlapSumRecFP.hpp"
#include "ThreeCenterOverlapSumRecFS.hpp"
#include "ThreeCenterOverlapSumRecDG.hpp"
#include "ThreeCenterOverlapSumRecDF.hpp"
#include "ThreeCenterOverlapSumRecDD.hpp"
#include "ThreeCenterOverlapSumRecDP.hpp"
#include "ThreeCenterOverlapSumRecDS.hpp"
#include "ThreeCenterOverlapSumRecPG.hpp"
#include "ThreeCenterOverlapSumRecPF.hpp"
#include "ThreeCenterOverlapSumRecPD.hpp"
#include "ThreeCenterOverlapSumRecPP.hpp"
#include "ThreeCenterOverlapSumRecPS.hpp"
#include "ThreeCenterOverlapSumRecSG.hpp"
#include "ThreeCenterOverlapSumRecSF.hpp"
#include "ThreeCenterOverlapSumRecSD.hpp"
#include "ThreeCenterOverlapSumRecSP.hpp"
#include "ThreeCenterOverlapSumRecSS.hpp"

namespace t3ovlfunc {

/// @brief Computes overlap integrals for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
compute(T&                               distributor,
        const CGtoBlock&                 bra_gto_block,
        const CGtoBlock&                 ket_gto_block,
        const std::pair<size_t, size_t>& bra_indices,
        const std::pair<size_t, size_t>& ket_indices,
        const bool                       bra_eq_ket) -> void
{
    const auto bra_angmom = bra_gto_block.angular_momentum();

    const auto ket_angmom = ket_gto_block.angular_momentum();

    if ((bra_angmom == 0) && (ket_angmom == 0))
    {
        t3ovlrec::comp_sum_overlap_ss(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_sp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        t3ovlrec::comp_sum_overlap_ps(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_pp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_sd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        t3ovlrec::comp_sum_overlap_ds(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_sf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        t3ovlrec::comp_sum_overlap_fs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_pd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_dp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_sg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        t3ovlrec::comp_sum_overlap_gs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_pf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_fp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_dd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_pg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_gp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_df(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_fd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_dg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_gd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_ff(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_fg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_gf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_gg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t3ovlfunc

#endif /* ThreeCenterOverlapFunc_hpp */
