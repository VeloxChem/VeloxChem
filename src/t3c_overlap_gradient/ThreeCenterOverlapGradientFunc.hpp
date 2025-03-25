#ifndef ThreeCenterOverlapGradientFunc_hpp
#define ThreeCenterOverlapGradientFunc_hpp

#include <array>

#include "GtoBlock.hpp"
#include "ThreeCenterOverlapGradientSumRecGG.hpp"
#include "ThreeCenterOverlapGradientSumRecGF.hpp"
#include "ThreeCenterOverlapGradientSumRecGD.hpp"
#include "ThreeCenterOverlapGradientSumRecGP.hpp"
#include "ThreeCenterOverlapGradientSumRecGS.hpp"
#include "ThreeCenterOverlapGradientSumRecFG.hpp"
#include "ThreeCenterOverlapGradientSumRecFF.hpp"
#include "ThreeCenterOverlapGradientSumRecFD.hpp"
#include "ThreeCenterOverlapGradientSumRecFP.hpp"
#include "ThreeCenterOverlapGradientSumRecFS.hpp"
#include "ThreeCenterOverlapGradientSumRecDG.hpp"
#include "ThreeCenterOverlapGradientSumRecDF.hpp"
#include "ThreeCenterOverlapGradientSumRecDD.hpp"
#include "ThreeCenterOverlapGradientSumRecDP.hpp"
#include "ThreeCenterOverlapGradientSumRecDS.hpp"
#include "ThreeCenterOverlapGradientSumRecPG.hpp"
#include "ThreeCenterOverlapGradientSumRecPF.hpp"
#include "ThreeCenterOverlapGradientSumRecPD.hpp"
#include "ThreeCenterOverlapGradientSumRecPP.hpp"
#include "ThreeCenterOverlapGradientSumRecPS.hpp"
#include "ThreeCenterOverlapGradientSumRecSG.hpp"
#include "ThreeCenterOverlapGradientSumRecSF.hpp"
#include "ThreeCenterOverlapGradientSumRecSD.hpp"
#include "ThreeCenterOverlapGradientSumRecSP.hpp"
#include "ThreeCenterOverlapGradientSumRecSS.hpp"

namespace g3ovlfunc {

/// @brief Computes overlap gradient integrals for given of pair basis functions blocks.
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
        g3ovlrec::comp_sum_overlap_gradient_ss(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        g3ovlrec::comp_sum_overlap_gradient_sp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        g3ovlrec::comp_sum_overlap_gradient_ps(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        g3ovlrec::comp_sum_overlap_gradient_pp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        g3ovlrec::comp_sum_overlap_gradient_sd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        g3ovlrec::comp_sum_overlap_gradient_ds(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        g3ovlrec::comp_sum_overlap_gradient_sf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        g3ovlrec::comp_sum_overlap_gradient_fs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        g3ovlrec::comp_sum_overlap_gradient_pd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        g3ovlrec::comp_sum_overlap_gradient_dp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        g3ovlrec::comp_sum_overlap_gradient_sg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        g3ovlrec::comp_sum_overlap_gradient_gs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        g3ovlrec::comp_sum_overlap_gradient_pf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        g3ovlrec::comp_sum_overlap_gradient_fp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        g3ovlrec::comp_sum_overlap_gradient_dd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        g3ovlrec::comp_sum_overlap_gradient_pg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        g3ovlrec::comp_sum_overlap_gradient_gp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        g3ovlrec::comp_sum_overlap_gradient_df(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        g3ovlrec::comp_sum_overlap_gradient_fd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        g3ovlrec::comp_sum_overlap_gradient_dg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        g3ovlrec::comp_sum_overlap_gradient_gd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        g3ovlrec::comp_sum_overlap_gradient_ff(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        g3ovlrec::comp_sum_overlap_gradient_fg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        g3ovlrec::comp_sum_overlap_gradient_gf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        g3ovlrec::comp_sum_overlap_gradient_gg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t3ovlfunc

#endif /* ThreeCenterOverlapGradientFunc_hpp */
